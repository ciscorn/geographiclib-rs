//! Local cartesian coordinates.

use crate::geocentric::Geocentric;
use crate::geomath;

/// A converter between geodetic coordinates and a local cartesian (East,
/// North, Up) coordinate system anchored at a chosen origin on the ellipsoid.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalCartesian {
    earth: Geocentric,
    lat0: f64,
    lon0: f64,
    h0: f64,
    x0: f64,
    y0: f64,
    z0: f64,
    r: [f64; 9],
}

impl LocalCartesian {
    /// Create a local cartesian converter with origin at `(lat0, lon0, h0)` on
    /// the WGS 84 ellipsoid.
    ///
    /// # Arguments
    ///
    /// * `lat0` - Latitude of the origin in degrees, in the range `[-90, 90]`.
    /// * `lon0` - Longitude of the origin in degrees.
    /// * `h0` - Height of the origin above the ellipsoid in meters.
    #[inline]
    pub fn wgs84(lat0: f64, lon0: f64, h0: f64) -> Self {
        Self::new(Geocentric::wgs84(), lat0, lon0, h0)
    }

    /// Create a local cartesian converter with origin at `(lat0, lon0, h0)` on
    /// the given ellipsoid.
    ///
    /// # Arguments
    ///
    /// * `earth` - Ellipsoid model.
    /// * `lat0` - Latitude of the origin in degrees, in the range `[-90, 90]`.
    /// * `lon0` - Longitude of the origin in degrees.
    /// * `h0` - Height of the origin above the ellipsoid in meters.
    #[inline]
    pub fn new(earth: Geocentric, lat0: f64, lon0: f64, h0: f64) -> Self {
        let mut lc = Self {
            earth,
            lat0: 0.0,
            lon0: 0.0,
            h0: 0.0,
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            r: [0.0; 9],
        };
        lc.reset(lat0, lon0, h0);
        lc
    }

    /// Reset the origin of the local cartesian frame to `(lat0, lon0, h0)`.
    ///
    /// # Arguments
    ///
    /// * `lat0` - Latitude of the origin in degrees, in the range `[-90, 90]`.
    /// * `lon0` - Longitude of the origin in degrees.
    /// * `h0` - Height of the origin above the ellipsoid in meters.
    #[inline]
    pub fn reset(&mut self, lat0: f64, lon0: f64, h0: f64) {
        let lat0 = geomath::lat_fix(lat0);
        let ((x0, y0, z0), r) = self.earth.forward_with_rotation(lat0, lon0, h0);
        self.lat0 = lat0;
        self.lon0 = lon0;
        self.h0 = h0;
        self.x0 = x0;
        self.y0 = y0;
        self.z0 = z0;
        self.r = r;
    }

    /// Latitude of the origin in degrees.
    #[inline]
    pub const fn latitude_origin(&self) -> f64 {
        self.lat0
    }

    /// Longitude of the origin in degrees.
    #[inline]
    pub const fn longitude_origin(&self) -> f64 {
        self.lon0
    }

    /// Height of the origin in meters above the ellipsoid.
    #[inline]
    pub const fn height_origin(&self) -> f64 {
        self.h0
    }

    /// Equatorial radius of the underlying ellipsoid in meters.
    #[inline]
    pub const fn equatorial_radius(&self) -> f64 {
        self.earth.equatorial_radius()
    }

    /// Flattening of the underlying ellipsoid.
    #[inline]
    pub const fn flattening(&self) -> f64 {
        self.earth.flattening()
    }

    /// Convert from geodetic coordinates to local cartesian coordinates.
    ///
    /// # Arguments
    ///
    /// * `lat` - Latitude in degrees, in the range `[-90, 90]`.
    /// * `lon` - Longitude in degrees.
    /// * `height` - Height above the ellipsoid in meters.
    ///
    /// # Returns
    ///
    /// `(x, y, z)` in meters in the local east, north, up frame anchored at
    /// the origin.
    #[inline]
    pub fn forward(&self, lat: f64, lon: f64, height: f64) -> (f64, f64, f64) {
        let (xc, yc, zc) = self.earth.forward(lat, lon, height);
        self.ecef_to_local(xc, yc, zc)
    }

    /// Convert from geodetic coordinates to local cartesian coordinates and
    /// return a rotation matrix.
    ///
    /// # Arguments
    ///
    /// * `lat` - Latitude in degrees, in the range `[-90, 90]`.
    /// * `lon` - Longitude in degrees.
    /// * `height` - Height above the ellipsoid in meters.
    ///
    /// # Returns
    ///
    /// A pair `((x, y, z), m)`, where `(x, y, z)` is the position in meters in
    /// the local east, north, up frame anchored at the origin, and `m` is a
    /// row-major 3x3 rotation matrix that maps a local ENU vector at
    /// `(lat, lon, height)` to the local cartesian frame.
    #[inline]
    pub fn forward_with_rotation(
        &self,
        lat: f64,
        lon: f64,
        height: f64,
    ) -> ((f64, f64, f64), [f64; 9]) {
        let ((xc, yc, zc), m) = self.earth.forward_with_rotation(lat, lon, height);
        let local = self.ecef_to_local(xc, yc, zc);
        (local, self.combine_rotation(&m))
    }

    /// Convert from local cartesian coordinates to geodetic coordinates.
    ///
    /// # Arguments
    ///
    /// * `x`, `y`, `z` - Local cartesian coordinates in meters in the east,
    ///   north, up frame anchored at the origin.
    ///
    /// # Returns
    ///
    /// `(lat, lon, height)`, where `lat` and `lon` are in degrees and `height`
    /// is in meters above the ellipsoid. If multiple geodetic solutions exist,
    /// the solution minimizing `abs(height)` is returned.
    #[inline]
    pub fn reverse(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let (xc, yc, zc) = self.local_to_ecef(x, y, z);
        self.earth.reverse(xc, yc, zc)
    }

    /// Convert from local cartesian coordinates to geodetic coordinates and
    /// return a rotation matrix.
    ///
    /// # Arguments
    ///
    /// * `x`, `y`, `z` - Local cartesian coordinates in meters in the east,
    ///   north, up frame anchored at the origin.
    ///
    /// # Returns
    ///
    /// A pair `((lat, lon, height), m)`, where `(lat, lon, height)` is the
    /// geodetic position (`lat` and `lon` in degrees, `height` in meters above
    /// the ellipsoid), and `m` is a row-major 3x3 rotation matrix that maps a
    /// local ENU vector at the returned geodetic position to the local
    /// cartesian frame.
    #[inline]
    pub fn reverse_with_rotation(&self, x: f64, y: f64, z: f64) -> ((f64, f64, f64), [f64; 9]) {
        let (xc, yc, zc) = self.local_to_ecef(x, y, z);
        let (lla, m) = self.earth.reverse_with_rotation(xc, yc, zc);
        (lla, self.combine_rotation(&m))
    }

    #[inline]
    fn ecef_to_local(&self, xc: f64, yc: f64, zc: f64) -> (f64, f64, f64) {
        let dx = xc - self.x0;
        let dy = yc - self.y0;
        let dz = zc - self.z0;
        let r = &self.r;
        // r maps ENU @ origin -> ECEF; r^T maps ECEF -> ENU @ origin.
        let x = r[0] * dx + r[3] * dy + r[6] * dz;
        let y = r[1] * dx + r[4] * dy + r[7] * dz;
        let z = r[2] * dx + r[5] * dy + r[8] * dz;
        (x, y, z)
    }

    #[inline]
    fn local_to_ecef(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let r = &self.r;
        let xc = r[0] * x + r[1] * y + r[2] * z + self.x0;
        let yc = r[3] * x + r[4] * y + r[5] * z + self.y0;
        let zc = r[6] * x + r[7] * y + r[8] * z + self.z0;
        (xc, yc, zc)
    }

    #[inline]
    fn combine_rotation(&self, m: &[f64; 9]) -> [f64; 9] {
        // Compute r^T * m, both 3x3 row-major. m maps ENU @ point -> ECEF, so
        // r^T * m maps ENU @ point -> ENU @ origin (the local cartesian frame).
        let r = &self.r;
        [
            r[0] * m[0] + r[3] * m[3] + r[6] * m[6],
            r[0] * m[1] + r[3] * m[4] + r[6] * m[7],
            r[0] * m[2] + r[3] * m[5] + r[6] * m[8],
            r[1] * m[0] + r[4] * m[3] + r[7] * m[6],
            r[1] * m[1] + r[4] * m[4] + r[7] * m[7],
            r[1] * m[2] + r[4] * m[5] + r[7] * m[8],
            r[2] * m[0] + r[5] * m[3] + r[8] * m[6],
            r[2] * m[1] + r[5] * m[4] + r[8] * m[7],
            r[2] * m[2] + r[5] * m[5] + r[8] * m[8],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn origin_maps_to_zero() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 100.0);
        let (x, y, z) = lc.forward(35.0, 139.0, 100.0);
        assert_abs_diff_eq!(x, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(y, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(z, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn reverse_at_zero_returns_origin() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 100.0);
        let (lat, lon, h) = lc.reverse(0.0, 0.0, 0.0);
        assert_abs_diff_eq!(lat, 35.0, epsilon = 1e-12);
        assert_abs_diff_eq!(lon, 139.0, epsilon = 1e-12);
        assert_abs_diff_eq!(h, 100.0, epsilon = 1e-7);
    }

    #[test]
    fn up_axis_increases_height() {
        // At the origin, moving along the local +z (up) axis by `dh` should
        // map back to the origin's lat/lon at height `h0 + dh`.
        let lc = LocalCartesian::wgs84(35.0, 139.0, 0.0);
        let (lat, lon, h) = lc.reverse(0.0, 0.0, 1234.5);
        assert_abs_diff_eq!(lat, 35.0, epsilon = 1e-12);
        assert_abs_diff_eq!(lon, 139.0, epsilon = 1e-12);
        assert_abs_diff_eq!(h, 1234.5, epsilon = 1e-7);
    }

    #[test]
    fn forward_at_equatorial_origin_axes() {
        // With origin at (0, 0, 0), local +x (east) aligns with ECEF +Y, +y
        // (north) aligns with ECEF +Z, and +z (up) aligns with ECEF +X.
        let lc = LocalCartesian::wgs84(0.0, 0.0, 0.0);

        let (x, y, z) = lc.forward(0.0, 1e-6, 0.0);
        assert!(x > 0.0 && y.abs() < 1e-3 && z.abs() < 1e-3);

        let (x, y, z) = lc.forward(1e-6, 0.0, 0.0);
        assert!(y > 0.0 && x.abs() < 1e-3 && z.abs() < 1e-3);

        let (x, y, z) = lc.forward(0.0, 0.0, 100.0);
        assert!(z > 0.0 && x.abs() < 1e-9 && y.abs() < 1e-9);
        assert_abs_diff_eq!(z, 100.0, epsilon = 1e-9);
    }

    #[test]
    fn roundtrip() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 50.0);
        for &(lat, lon, h) in &[
            (35.5, 139.5, 200.0),
            (34.0, 140.5, -10.0),
            (-10.0, -170.0, 0.0),
            (89.0, 0.0, 1000.0),
        ] {
            let (x, y, z) = lc.forward(lat, lon, h);
            let (lat2, lon2, h2) = lc.reverse(x, y, z);
            assert_abs_diff_eq!(lat, lat2, epsilon = 1e-10);
            assert_abs_diff_eq!(lon, lon2, epsilon = 1e-10);
            assert_abs_diff_eq!(h, h2, epsilon = 1e-7);
        }
    }

    #[test]
    fn forward_with_rotation_matches_forward() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 0.0);
        let (lat, lon, h) = (35.5, 139.5, 200.0);
        let plain = lc.forward(lat, lon, h);
        let (with_rot, _) = lc.forward_with_rotation(lat, lon, h);
        assert_abs_diff_eq!(plain.0, with_rot.0, epsilon = 1e-12);
        assert_abs_diff_eq!(plain.1, with_rot.1, epsilon = 1e-12);
        assert_abs_diff_eq!(plain.2, with_rot.2, epsilon = 1e-12);
    }

    #[test]
    fn rotation_at_origin_is_identity() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 100.0);
        let (_, m) = lc.forward_with_rotation(35.0, 139.0, 100.0);
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        for (a, b) in m.iter().zip(identity.iter()) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-12);
        }
    }

    #[test]
    fn rotation_maps_local_enu_to_local_cartesian() {
        // Pick a point that is not the origin, get its local ENU -> local
        // cartesian rotation, and verify that the unit "up" vector at the
        // queried point points away from the earth's center as seen in the
        // local cartesian frame (i.e. has a positive radial component).
        let lc = LocalCartesian::wgs84(35.0, 139.0, 0.0);
        let (pos, m) = lc.forward_with_rotation(35.5, 139.5, 0.0);
        // The "up" vector (0, 0, 1) at the queried point in local cartesian.
        let up_local = (m[2], m[5], m[8]);
        // Should be roughly parallel to the position vector from the origin
        // (since the origin is on the ellipsoid and so is the queried point).
        let pos_norm = (pos.0 * pos.0 + pos.1 * pos.1 + pos.2 * pos.2).sqrt();
        let dot = (up_local.0 * pos.0 + up_local.1 * pos.1 + up_local.2 * pos.2) / pos_norm;
        // The two vectors should have a substantial positive component in
        // common (the queried point is "outward" relative to the origin).
        assert!(dot > 0.0);
    }

    #[test]
    fn reverse_with_rotation_inverts_forward_with_rotation() {
        let lc = LocalCartesian::wgs84(35.0, 139.0, 50.0);
        let (lat, lon, h) = (35.5, 139.5, 200.0);
        let ((x, y, z), forward_rot) = lc.forward_with_rotation(lat, lon, h);
        let ((lat2, lon2, h2), reverse_rot) = lc.reverse_with_rotation(x, y, z);
        assert_abs_diff_eq!(lat, lat2, epsilon = 1e-10);
        assert_abs_diff_eq!(lon, lon2, epsilon = 1e-10);
        assert_abs_diff_eq!(h, h2, epsilon = 1e-7);
        for (f, r) in forward_rot.iter().zip(reverse_rot.iter()) {
            assert_abs_diff_eq!(f, r, epsilon = 1e-12);
        }
    }

    #[test]
    fn reset_changes_origin() {
        let mut lc = LocalCartesian::wgs84(0.0, 0.0, 0.0);
        let (x_before, _, _) = lc.forward(0.0, 1e-6, 0.0);
        lc.reset(35.0, 139.0, 0.0);
        assert_eq!(lc.latitude_origin(), 35.0);
        assert_eq!(lc.longitude_origin(), 139.0);
        assert_eq!(lc.height_origin(), 0.0);
        let (x_after, y_after, z_after) = lc.forward(35.0, 139.0, 0.0);
        assert_abs_diff_eq!(x_after, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(y_after, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(z_after, 0.0, epsilon = 1e-9);
        // Sanity: the pre-reset and post-reset frames are different.
        assert_ne!(x_before, x_after);
    }

    #[test]
    fn out_of_range_origin_latitude_propagates_nan() {
        let lc = LocalCartesian::wgs84(91.0, 0.0, 0.0);
        assert!(lc.latitude_origin().is_nan());
        let (x, y, z) = lc.forward(0.0, 0.0, 0.0);
        assert!(x.is_nan());
        assert!(y.is_nan());
        assert!(z.is_nan());
    }

    #[test]
    fn custom_ellipsoid() {
        let earth = Geocentric::new(6_400_000.0, 1.0 / 300.0);
        let lc = LocalCartesian::new(earth, 35.0, 139.0, 0.0);
        assert_eq!(lc.equatorial_radius(), 6_400_000.0);
        assert_eq!(lc.flattening(), 1.0 / 300.0);
        let (x, y, z) = lc.forward(35.0, 139.0, 0.0);
        assert_abs_diff_eq!(x, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(y, 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(z, 0.0, epsilon = 1e-9);
    }
}
