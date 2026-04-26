//! Geocentric calculations.
//!
//! Geocentric coordinates are also known as earth-centered, earth-fixed (ECEF)
//! coordinates. The origin is at the center of the earth. The Z axis passes
//! through the north pole, the X axis passes through latitude 0 and longitude
//! 0, and the Y axis completes a right-handed coordinate system.

use core::f64::consts::FRAC_PI_6;
use std::sync;

use crate::geodesic::{WGS84_A, WGS84_F};
use crate::geomath;

/// A converter between geodetic coordinates and geocentric coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Geocentric {
    a: f64,
    f: f64,
    e_sq: f64,
    e_quad: f64,
    max_rad: f64,
    one_minus_e_sq: f64,
}

static WGS84_GEOCENTRIC: sync::OnceLock<Geocentric> = sync::OnceLock::new();

impl Geocentric {
    /// Create a geocentric converter for the WGS 84 ellipsoid.
    #[inline]
    pub fn wgs84() -> Self {
        *WGS84_GEOCENTRIC.get_or_init(|| Self::new(WGS84_A, WGS84_F))
    }

    /// Create a geocentric converter for an ellipsoid of revolution.
    ///
    /// `a` is the equatorial radius in meters. `f` is the flattening of the
    /// ellipsoid. Setting `f = 0` gives a sphere; negative `f` gives a prolate
    /// ellipsoid.
    ///
    /// # Panics
    /// Panics if `a` is not positive and finite, or if `(1 - f) * a` is not
    /// positive and finite.
    #[inline]
    pub fn new(a: f64, f: f64) -> Self {
        assert!(a.is_finite() && a > 0.0, "semi-major axis must be positive");
        assert!(
            f.is_finite() && ((1. - f) * a).is_finite() && (1. - f) * a > 0.,
            "polar semi-axis must be positive"
        );
        Self::new_unchecked(a, f)
    }

    #[inline]
    const fn new_unchecked(a: f64, f: f64) -> Self {
        let e_sq = f * (2.0 - f);
        let one_minus_f = 1.0 - f;
        Self {
            a,
            f,
            e_sq,
            e_quad: e_sq * e_sq,
            max_rad: 2.0 * a / f64::EPSILON,
            one_minus_e_sq: one_minus_f * one_minus_f,
        }
    }

    /// Equatorial radius in meters.
    #[inline]
    pub const fn equatorial_radius(&self) -> f64 {
        self.a
    }

    /// Flattening of the ellipsoid.
    #[inline]
    pub const fn flattening(&self) -> f64 {
        self.f
    }

    /// Convert from geodetic coordinates to geocentric coordinates.
    ///
    /// `lat` and `lon` are in degrees, and `height` is in meters above the
    /// ellipsoid. `lat` should be in the range `[-90, 90]`.
    ///
    /// Returns `(x, y, z)` in meters.
    #[inline]
    pub fn forward(&self, lat: f64, lon: f64, height: f64) -> (f64, f64, f64) {
        let (sin_phi, cos_phi) = geomath::sincosd(geomath::lat_fix(lat));
        let (sin_lam, cos_lam) = geomath::sincosd(lon);
        let n = if self.e_sq == 0.0 {
            self.a
        } else {
            self.a / (1. - self.e_sq * sin_phi * sin_phi).sqrt()
        };
        let r = (n + height) * cos_phi;
        let x = r * cos_lam;
        let y = r * sin_lam;
        let z = (n * (1. - self.e_sq) + height) * sin_phi;
        (x, y, z)
    }

    /// Convert from geodetic coordinates to geocentric coordinates and return a
    /// rotation matrix.
    ///
    /// The returned matrix is row-major and maps a local east, north, up (ENU)
    /// vector at `(lat, lon, height)` to a geocentric ECEF vector.
    #[inline]
    pub fn forward_with_rotation(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
    ) -> ((f64, f64, f64), [f64; 9]) {
        let xyz = self.forward(lat, lon, height);
        let (sin_lam, cos_lam) = geomath::sincosd(lon);
        let (sin_phi, cos_phi) = geomath::sincosd(geomath::lat_fix(lat));
        let rotation = rotation_matrix(sin_phi, cos_phi, sin_lam, cos_lam);
        (xyz, rotation)
    }

    /// Convert from geocentric coordinates to geodetic coordinates.
    ///
    /// `x`, `y`, and `z` are geocentric coordinates in meters.
    ///
    /// Returns `(lat, lon, height)`, where `lat` and `lon` are in degrees and
    /// `height` is in meters above the ellipsoid. If multiple geodetic
    /// solutions exist, the solution minimizing `abs(height)` is returned.
    ///
    /// This method is based on Vermeille's method with GeographicLib's
    /// robustness improvement.
    #[inline]
    pub fn reverse(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let rho = x.hypot(y);
        let height_to_center = rho.hypot(z);

        if height_to_center > self.max_rad {
            // Match GeographicLib's overflow guard: for extremely distant points,
            // treating the earth as a point gives an adequate height approximation.
            let half_x = x / 2.0;
            let half_y = y / 2.0;
            let half_z = z / 2.0;
            let half_rho = half_x.hypot(half_y);
            let (sin_lam, cos_lam) = if half_rho != 0.0 {
                (half_y / half_rho, half_x / half_rho)
            } else {
                (0.0, 1.0)
            };
            let half_h = half_z.hypot(half_rho);
            let sin_phi = half_z / half_h;
            let cos_phi = half_rho / half_h;
            return (
                geomath::atan2d(sin_phi, cos_phi),
                geomath::atan2d(sin_lam, cos_lam),
                height_to_center,
            );
        }

        if self.e_sq == 0.0 {
            let phi = if height_to_center == 0.0 {
                90.0
            } else {
                geomath::atan2d(z, rho)
            };
            let lam = geomath::atan2d(y, x);
            return (phi, lam, height_to_center - self.a);
        }

        let prolate = self.e_sq < 0.;
        let e_abs = self.e_sq.abs();
        let pa = rho / self.a;
        let za = z / self.a;
        let mut p = pa * pa;
        let mut q = self.one_minus_e_sq * za * za;
        let r = (p + q - self.e_quad) / 6.;
        if prolate {
            core::mem::swap(&mut p, &mut q);
        }
        let r3 = r * r * r;
        let e_quad_p_q = self.e_quad * p * q;

        let evol = 8. * r3 + e_quad_p_q;

        let (sin_phi, cos_phi, h) = if evol > 0. || q != 0. {
            let u = if evol > 0. {
                let l = (evol.sqrt() + e_quad_p_q.sqrt()).cbrt();
                (3. * r * r) / (2. * l * l) + 0.5 * (l + r / l) * (l + r / l)
            } else {
                let t = 2. / 3. * e_quad_p_q.sqrt().atan2((-evol).sqrt() + (-8. * r3).sqrt());
                -4. * r * (t).sin() * (FRAC_PI_6 + t).cos()
            };
            let v = (u * u + self.e_quad * q).sqrt();
            let uv = if u < 0. {
                self.e_quad * q / (v - u)
            } else {
                u + v
            };
            let w = (e_abs * (uv - q) / (2. * v)).max(0.);
            let k = uv / ((w * w + uv).sqrt() + w);
            let (k1, k2) = if prolate {
                (k - self.e_sq, k)
            } else {
                (k, k + self.e_sq)
            };
            let d = k1 * rho / k2;
            let h = (1. - self.one_minus_e_sq / k1) * d.hypot(z);
            let big_h = (z / k1).hypot(rho / k2);
            let sin_phi = (z / k1) / big_h;
            let cos_phi = (rho / k2) / big_h;
            (sin_phi, cos_phi, h)
        } else {
            let zz = ((if prolate { p } else { self.e_quad - p }) / self.one_minus_e_sq).sqrt();
            let xx = (if prolate { self.e_quad - p } else { p }).sqrt();
            let big_h = zz.hypot(xx);
            let mut sin_phi = zz / big_h;
            let cos_phi = xx / big_h;
            if z < 0. {
                sin_phi = -sin_phi;
            }
            let h = -self.a * (if prolate { 1. } else { self.one_minus_e_sq }) * big_h / e_abs;
            (sin_phi, cos_phi, h)
        };

        let lam = if rho != 0.0 {
            geomath::atan2d(y / rho, x / rho)
        } else {
            0.0
        };

        (geomath::atan2d(sin_phi, cos_phi), lam, h)
    }

    /// Convert from geocentric coordinates to geodetic coordinates and return a
    /// rotation matrix.
    ///
    /// The returned matrix is row-major and maps a local east, north, up (ENU)
    /// vector at the returned geodetic position to a geocentric ECEF vector.
    #[inline]
    pub fn reverse_with_rotation(&self, x: f64, y: f64, z: f64) -> ((f64, f64, f64), [f64; 9]) {
        let lla = self.reverse(x, y, z);
        let (sin_lam, cos_lam) = geomath::sincosd(lla.1);
        let (sin_phi, cos_phi) = geomath::sincosd(lla.0);
        let rotation = rotation_matrix(sin_phi, cos_phi, sin_lam, cos_lam);
        (lla, rotation)
    }
}

#[inline]
fn rotation_matrix(sin_phi: f64, cos_phi: f64, sin_lam: f64, cos_lam: f64) -> [f64; 9] {
    [
        -sin_lam,
        -cos_lam * sin_phi,
        cos_lam * cos_phi,
        cos_lam,
        -sin_lam * sin_phi,
        sin_lam * cos_phi,
        0.,
        cos_phi,
        sin_phi,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use std::io::BufRead;

    const GEOCENTRIC_WGS84_CPP_TEST_PATH: &str = "test_fixtures/geocentric_wgs84_cpp.dat";

    #[test]
    fn roundtrip() {
        let a = 6378137.;
        let inv_f = 298.257223563;
        let f = 1. / inv_f;
        let b = a * (1. - f);
        let earth = Geocentric::new(a, f);

        {
            let (lat, lon, height) = (37., 140., 50.);
            let (x, y, z) = earth.forward(lat, lon, height);
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((lon - lon2).abs() < 1e-10);
            assert!((height - height2).abs() < 1e-7);
        }

        {
            let (lat, lon, height) = (74.58501644931525, 45., -6344866.234164982);
            let (x, y, z) = earth.forward(lat, lon, height);
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((lon - lon2).abs() < 1e-10);
            assert!((height - height2).abs() < 1e-7);
        }

        {
            let (lat, lon, height) = (88.10828645, 120., -6356728.972246517);
            let (x, y, _) = earth.forward(lat, lon, height);
            let z = 0.;
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((lon - lon2).abs() < 1e-10);
            assert!((height - height2).abs() < 30.);
        }

        {
            let (lat, lon, height) = (-88.10828645, 120., -6356728.972246517);
            let (x, y, _) = earth.forward(lat, lon, height);
            let z = -f64::MIN_POSITIVE;
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((lon - lon2).abs() < 1e-10);
            assert!((height - height2).abs() < 30.);
        }

        {
            let (lat, lon, height) = (90., 0., b);
            let (x, y, z) = earth.forward(lat, lon, height);
            let (lat2, lon2, height2) = earth.reverse(x, y, z);
            assert!((lat - lat2).abs() < 1e-9);
            assert!((lon - lon2).abs() < 1e-9);
            assert!((height - height2).abs() < 1e-7);
        }
    }

    #[test]
    fn to_geocentric() {
        let a = 6378137.;
        let inv_f = 298.257223563;
        let f = 1. / inv_f;
        let b = a * (1. - f);
        let earth = Geocentric::new(a, f);

        {
            let (x, y, z) = earth.forward(37., 140., 50.);
            assert!((x - -3906851.9770472576).abs() < 1e-9);
            assert!((y - 3278238.0530045824).abs() < 1e-9);
            assert!((z - 3817423.251099322).abs() < 1e-9);
        }

        {
            let height = 150.;
            let (x, y, z) = earth.forward(90., 123., height);
            assert!((x - 0.).abs() < 1e-9);
            assert!((y - 0.).abs() < 1e-9);
            assert!((z - (b + height)).abs() < 1e-9);
        }

        {
            let height = 100.;
            let (x, y, z) = earth.forward(0., 0., height);
            assert!((x - (a + height)).abs() < 1e-9);
            assert!((y - 0.).abs() < 1e-9);
            assert!((z - 0.).abs() < 1e-9);
        }
    }

    #[test]
    fn sphere_reverse() {
        let earth = Geocentric::new(10., 0.);
        assert_eq!(earth.flattening(), 0.);

        {
            let (lat, lon, height) = earth.reverse(0., 0., 0.);
            assert_eq!(lat, 90.);
            assert_eq!(lon, 0.);
            assert_eq!(height, -10.);
        }

        {
            let (lat, lon, height) = earth.reverse(10., 0., 0.);
            assert_eq!(lat, 0.);
            assert_eq!(lon, 0.);
            assert_eq!(height, 0.);
        }
    }

    #[test]
    fn sphere_forward() {
        let earth = Geocentric::new(10., 0.);
        let (x, y, z) = earth.forward(0., 0., 5.);

        assert_eq!(x, 15.);
        assert_eq!(y, 0.);
        assert_eq!(z, 0.);
    }

    #[test]
    fn forward_out_of_range_latitude_returns_nan() {
        let earth = Geocentric::wgs84();

        let (x, y, z) = earth.forward(91., 0., 0.);
        assert!(x.is_nan());
        assert!(y.is_nan());
        assert!(z.is_nan());

        let ((x, y, z), rotation) = earth.forward_with_rotation(-91., 0., 0.);
        assert!(x.is_nan());
        assert!(y.is_nan());
        assert!(z.is_nan());
        assert!(rotation.iter().any(|value| value.is_nan()));
    }

    #[test]
    fn geographiclib_cartconvert_reverse_cases() {
        {
            let earth = Geocentric::new(6.4e6, 1. / 100.);
            let (lat, lon, height) = earth.reverse(10e3, 0., 1e3);
            assert!((lat - 85.57).abs() < 0.01);
            assert!((lon - 0.).abs() < 1e-12);
            assert!((height - -6334614.).abs() < 1.);
        }

        {
            let earth = Geocentric::new(6.4e6, -1. / 100.);
            let (lat, lon, height) = earth.reverse(1e3, 0., 10e3);
            assert!((lat - 4.42).abs() < 0.01);
            assert!((lon - 0.).abs() < 1e-12);
            assert!((height - -6398614.).abs() < 1.);
        }
    }

    #[test]
    fn reverse_far_point_uses_geographiclib_overflow_guard() {
        let earth = Geocentric::wgs84();
        let huge = 1.0e23_f64;

        let (lat, lon, height) = earth.reverse(huge, huge, huge);
        let expected_height = huge.hypot(huge).hypot(huge);

        assert!((lat - 35.264389682754654).abs() < 1e-12);
        assert!((lon - 45.).abs() < 1e-12);
        assert!((height - expected_height).abs() / expected_height < 1e-15);

        let (lat, lon, height) = earth.reverse(f64::MAX, f64::MAX, f64::MAX);

        assert!((lat - 35.264389682754654).abs() < 1e-12);
        assert!((lon - 45.).abs() < 1e-12);
        assert!(height.is_infinite());

        let (lat, lon, height) = earth.reverse(0., 0., huge);

        assert_eq!(lat, 90.);
        assert_eq!(lon, 0.);
        assert_eq!(height, huge);
    }

    #[test]
    fn prolate_roundtrip() {
        let earth = Geocentric::new(6.4e6, -1. / 100.);
        let (lat, lon, height) = (35., 123., 1000.);

        let (x, y, z) = earth.forward(lat, lon, height);
        let (lat2, lon2, height2) = earth.reverse(x, y, z);

        assert!((lat - lat2).abs() < 1e-10);
        assert!((lon - lon2).abs() < 1e-10);
        assert!((height - height2).abs() < 1e-7);
    }

    #[test]
    fn forward_with_rotation() {
        let earth = Geocentric::wgs84();
        let ((x, y, z), m) = earth.forward_with_rotation(0., 0., 0.);

        assert!((x - earth.equatorial_radius()).abs() < 1e-9);
        assert_eq!(y, 0.);
        assert_eq!(z, 0.);
        assert_eq!(m, [-0., -0., 1., 1., -0., 0., 0., 1., 0.]);
    }

    #[test]
    fn reverse_with_rotation_matches_forward_rotation() {
        let earth = Geocentric::wgs84();
        let (lat, lon, height) = (37., 140., 50.);
        let ((x, y, z), forward_rotation) = earth.forward_with_rotation(lat, lon, height);
        let ((lat2, lon2, height2), reverse_rotation) = earth.reverse_with_rotation(x, y, z);

        assert!((lat - lat2).abs() < 1e-10);
        assert!((lon - lon2).abs() < 1e-10);
        assert!((height - height2).abs() < 1e-7);
        for (forward, reverse) in forward_rotation.iter().zip(reverse_rotation.iter()) {
            assert!((forward - reverse).abs() < 1e-12);
        }
    }

    #[test]
    fn geocentric_cpp_compatibility() {
        let earth = Geocentric::wgs84();
        // Generated from C++ GeographicLib Geocentric with std::setprecision(17).
        let file = std::fs::File::open(GEOCENTRIC_WGS84_CPP_TEST_PATH)
            .expect("failed to open geocentric_wgs84_cpp.dat");
        let reader = std::io::BufReader::new(file);

        for (i, line) in reader.lines().enumerate() {
            let line_num = i + 1;
            let line = line.expect("failed to read GeocentricTest line");
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|value| {
                    value
                        .parse::<f64>()
                        .expect("failed to parse GeocentricTest value")
                })
                .collect();
            assert_eq!(values.len(), 9, "line {line_num}");

            let (x, y, z) = earth.forward(values[0], values[1], values[2]);
            assert_relative_eq!(x, values[3], epsilon = 1e-8, max_relative = 1e-14);
            assert_relative_eq!(y, values[4], epsilon = 1e-8, max_relative = 1e-14);
            assert_relative_eq!(z, values[5], epsilon = 1e-8, max_relative = 1e-14);

            let (lat, lon, height) = earth.reverse(values[3], values[4], values[5]);
            assert_abs_diff_eq!(lat, values[6], epsilon = 1e-12);
            assert_abs_diff_eq!(lon, values[7], epsilon = 1e-12);
            assert_relative_eq!(height, values[8], epsilon = 1e-7, max_relative = 1e-14);
        }
    }
}
