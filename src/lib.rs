use chrono::prelude::*;
use chrono::DateTime;

use thiserror::Error;

mod sgp4_sys;

#[derive(Debug, Error)]
pub enum Error {
    #[error("TLE was malformed: {0}")]
    MalformedTwoLineElement(String),
    #[error("Error in SGP4 propagator")]
    PropagationError,
    #[error("{0}")]
    UnknownError(String),
}

type Result<T> = std::result::Result<T, Error>;

/// A state vector, in the TEME-ECI coordinate frame, for an orbiting body.
///
/// To obtain the state of an object at a specific time, use the propagation functions provided by
/// [TwoLineElement].
#[derive(Debug)]
pub struct StateVector {
    pub epoch: DateTime<Utc>,

    /// The satellite position in km.
    pub position: [f64; 3],

    /// The satellite velocity in km/s.
    pub velocity: [f64; 3],
}

const TLE_LINE_LENGTH: usize = 69;

pub struct TwoLineElement {
    elements: sgp4_sys::OrbitalElementSet,
}

impl TwoLineElement {
    /// Create a validated TwoLineElement from a string.
    pub fn new(line1: &str, line2: &str) -> Result<TwoLineElement> {
        let line1 = line1.trim();
        let line2 = line2.trim();

        if line1.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement(format!(
                "Line 1 is the wrong length. Expected {}, but got {}",
                TLE_LINE_LENGTH,
                line1.len()
            )));
        }

        if line2.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement(format!(
                "Line 2 is the wrong length. Expected {}, but got {}",
                TLE_LINE_LENGTH,
                line2.len()
            )));
        }

        let elements = sgp4_sys::to_orbital_elements(
            line1,
            line2,
            sgp4_sys::RunType::Verification,
            sgp4_sys::OperationMode::Improved,
            sgp4_sys::GravitationalConstant::Wgs84,
        )
        .map_err(|e| Error::MalformedTwoLineElement(format!("{:?}", e)))?;

        Ok(TwoLineElement { elements })
    }

    /// Create a TwoLineElement from a string containing both lines, and optionally a header line.
    pub fn from_lines(combined_lines: &str) -> Result<TwoLineElement> {
        let lines: Vec<_> = {
            let mut ls: Vec<_> = combined_lines
                .split('\n')
                .filter(|s| !s.is_empty())
                .collect();
            if ls.len() == 3 {
                ls.split_off(1)
            } else if ls.len() == 2 {
                ls
            } else {
                return Err(Error::MalformedTwoLineElement(format!(
                    "Expected two lines, got {}",
                    ls.len()
                )));
            }
        };
        TwoLineElement::new(&lines[0], &lines[1])
    }

    /// Get the epoch of a TwoLineElement.
    pub fn epoch(&self) -> Result<DateTime<Utc>> {
        Ok(self.elements.epoch())
    }

    /// Propagate a TwoLineElement to the given time to obtain a state vector for the object.
    pub fn propagate_to(&self, t: DateTime<Utc>) -> Result<StateVector> {
        let tle_epoch = self.elements.epoch();
        // TODO: determine correct behaviour for negative prop
        // assert!(t >= tle_epoch);

        let min_since_epoch = (t - tle_epoch).num_milliseconds() as f64 / 60_000.;

        let (r, v) = sgp4_sys::run_sgp4(
            self.elements,
            sgp4_sys::GravitationalConstant::Wgs84,
            min_since_epoch,
        )
        .map_err(|_e| Error::PropagationError)?;

        Ok(StateVector {
            epoch: t,
            position: r.to_owned(),
            velocity: v.to_owned(),
        })
    }
}

/// Wrapper type representing a Julian day.
///
/// This is the number of days since the start of the Julian astronomical calendar in 4713 BC, used
/// to provide a consistent time reference for astronomical calculations.
pub struct JulianDay(f64);

impl From<DateTime<Utc>> for JulianDay {
    fn from(d: DateTime<Utc>) -> Self {
        JulianDay(sgp4_sys::datetime_to_julian_day(d) as f64)
    }
}

impl From<JulianDay> for DateTime<Utc> {
    fn from(jd: JulianDay) -> Self {
        sgp4_sys::julian_day_to_datetime(jd.0)
    }
}

/// Wrapper type representing the angular form of Greenwich Mean Sidereal Time.
///
/// This is primarily used to account for the Earth's rotation during conversion between fixed and
/// inertial coordinate frames. Note that this is an angle measured in radians, and not a "time" as
/// such. The value may range from 0 to 2π.
pub struct GreenwichMeanSiderealTime(f64);

impl GreenwichMeanSiderealTime {
    pub fn as_radians(&self) -> f64 {
        self.0
    }
}

impl From<DateTime<Utc>> for GreenwichMeanSiderealTime {
    fn from(d: DateTime<Utc>) -> Self {
        GreenwichMeanSiderealTime(sgp4_sys::datetime_to_gstime(d) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use chrono::Duration;
    use float_cmp::approx_eq;

    fn vecs_eq(l: &[f64; 3], r: &[f64; 3]) -> bool {
        approx_eq!(f64, l[0], r[0])
        && approx_eq!(f64, l[1], r[1])
        && approx_eq!(f64, l[2], r[2])
    }

    #[test]
    fn test_simple_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch)?;
        let s2 = tle.propagate_to(epoch + Duration::hours(1))?;

        assert!(!vecs_eq(&s1.position, &s2.position));
        assert!(!vecs_eq(&s1.velocity, &s2.velocity));

        Ok(())
    }

    #[test]
    fn test_negative_time_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch)?;
        let s2 = tle.propagate_to(epoch - Duration::days(30))?;


        assert!(!vecs_eq(&s1.position, &s2.position));
        assert!(!vecs_eq(&s1.velocity, &s2.velocity));

        Ok(())
    }

    #[test]
    fn test_tle_from_lines() -> Result<()> {
        let lines = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let _tle = TwoLineElement::from_lines(lines)?;
        Ok(())
    }

    #[test]
    fn test_tle_from_lines_with_header() -> Result<()> {
        let lines = "ISS (ZARYA)
                     1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let _tle = TwoLineElement::from_lines(lines)?;
        Ok(())
    }

    #[test]
    fn test_tle_from_lines_with_surrounding_whitespace() -> Result<()> {
        let lines = "\nISS (ZARYA)
                     1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
                     2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767\n";

        let _tle = TwoLineElement::from_lines(lines)?;

        Ok(())
    }

    #[test]
    fn test_julian_day_identity() {
        let t = Utc.ymd(2020, 1, 1).and_hms(0, 0, 0);
        assert_eq!(DateTime::<Utc>::from(JulianDay::from(t)), t);
    }

    #[test]
    fn test_gmst_conversion() {
        let t = Utc.ymd(2020, 1, 1).and_hms(0, 0, 0);
        let a: f64 = 100.1218209532; // GMST for 2020-01-01T00:00:00 in degrees
        let a_rad = a.to_radians();
        assert!(sgp4_sys::close(
            GreenwichMeanSiderealTime::from(t).as_radians(),
            a_rad
        ));
    }
}
