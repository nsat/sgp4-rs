//! # SGP-4 System
//!
//! This module contains FFI bindings to the Vallado C++ implementation of SGP-4. These bindings
//! rely heavily on `unsafe` code, and it is recommended to instead use the other, high level
//! bindings provided by the `sgp4` crate.

use std::ffi::{CString, NulError};
use std::os::raw::{c_char, c_double, c_int, c_long};

use chrono::prelude::*;

use thiserror::Error;

#[derive(Debug, Error, PartialEq)]
pub enum Error {
    #[error(transparent)]
    CStringNul(#[from] NulError),
    #[error("Unknown failure in SGP4 propagator")]
    Unknown,
    #[error("Eccentricity out of bounds for mean elements")]
    InvalidEccentricity,
    #[error("Mean motion must be positive")]
    NegativeMeanMotion,
    #[error("Eccentricity out of bounds for pert elements")]
    PertEccentricity,
    #[error("Semi-latus rectum must be positive")]
    NegativeSemiLatus,
    #[error("Epoch elements are sub-orbital")]
    SubOrbital,
    #[error("Satellite has decayed")]
    SatelliteDecay,
}

#[allow(dead_code)]
#[derive(Default)]
pub(crate) enum RunType {
    #[default]
    Verification,
    Catalog,
    Manual(InputType),
}

impl RunType {
    fn to_char(&self) -> c_char {
        use RunType::*;

        match self {
            Verification => 'v' as c_char,
            Catalog => 'c' as c_char,
            Manual(_) => 'm' as c_char,
        }
    }

    fn to_input_char(&self) -> c_char {
        use RunType::*;

        match self {
            Manual(it) => it.to_char(),
            _ => 'v' as c_char, // Value is ignored by SGP-4. See sgp4io.cpp:196.
        }
    }
}

#[allow(dead_code)]
pub(crate) enum InputType {
    Epoch,
    MinutesFromEpoch,
    DayOfYear,
}

impl InputType {
    fn to_char(&self) -> c_char {
        use InputType::*;
        match self {
            Epoch => 'e' as c_char,
            MinutesFromEpoch => 'm' as c_char,
            DayOfYear => 'd' as c_char,
        }
    }
}

#[allow(dead_code)]
#[derive(Default)]
pub(crate) enum OperationMode {
    AirForceSpaceCenter,
    #[default]
    Improved,
}

impl OperationMode {
    fn to_char(&self) -> c_char {
        use OperationMode::*;

        match self {
            AirForceSpaceCenter => 'a' as c_char,
            Improved => 'i' as c_char,
        }
    }
}

const EPSILON: f64 = 0.000_000_1;

#[repr(C)]
#[allow(dead_code, non_camel_case_types)]
#[derive(Copy, Clone, Debug)]
pub(crate) enum GravitationalConstant {
    Wgs72Old,
    Wgs72,
    Wgs84,
}

#[repr(C)]
#[allow(dead_code)]
#[derive(Default, Clone, Copy, Debug)]
pub(crate) struct OrbitalElementSet {
    catalog_number: c_long, // satnum
    epoch_year: c_int,      // epochyr

    // Unused?
    epochtynumrev: c_int,

    error: c_int,

    operation_mode: c_char,
    init: c_char,
    method: c_char,

    /* Near Earth */
    isimp: c_int,
    aycof: c_double,
    con41: c_double,
    cc1: c_double,
    cc4: c_double,
    cc5: c_double,
    d2: c_double,
    d3: c_double,
    d4: c_double,

    delmo: c_double,
    eta: c_double,
    argpdot: c_double,
    omgcof: c_double,
    sinmao: c_double,
    t: c_double,
    t2cof: c_double,
    t3cof: c_double,

    t4cof: c_double,
    t5cof: c_double,
    x1mth2: c_double,
    x7thm1: c_double,
    mdot: c_double,
    nodedot: c_double,
    xlcof: c_double,
    xmcof: c_double,

    nodecf: c_double,

    /* Deep Space */
    irez: c_int,
    d2201: c_double,
    d2211: c_double,
    d3210: c_double,
    d3222: c_double,
    d4410: c_double,
    d4422: c_double,
    d5220: c_double,
    d5232: c_double,

    d5421: c_double,
    d5433: c_double,
    dedt: c_double,
    del1: c_double,
    del2: c_double,
    del3: c_double,
    didt: c_double,
    dmdt: c_double,

    dnodt: c_double,
    domdt: c_double,
    e3: c_double,
    ee2: c_double,
    peo: c_double,
    pgho: c_double,
    pho: c_double,
    pinco: c_double,

    plo: c_double,
    se2: c_double,
    se3: c_double,
    sgh2: c_double,
    sgh3: c_double,
    sgh4: c_double,
    sh2: c_double,
    sh3: c_double,

    si2: c_double,
    si3: c_double,
    sl2: c_double,
    sl3: c_double,
    sl4: c_double,
    gsto: c_double,
    xfact: c_double,
    xgh2: c_double,

    xgh3: c_double,
    xgh4: c_double,
    xh2: c_double,
    xh3: c_double,
    xi2: c_double,
    xi3: c_double,
    xl2: c_double,
    xl3: c_double,

    xl4: c_double,
    xlamo: c_double,
    zmol: c_double,
    zmos: c_double,
    atime: c_double,
    xli: c_double,
    xni: c_double,

    // The following fields are used in TLEs and provide the "classical" orbital elements.
    pub(crate) semi_major_axis: c_double,       // a
    pub(crate) altitude_of_periapsis: c_double, // altp
    pub(crate) altitude_of_apoapsis: c_double,  // alta

    epoch_days: c_double,                    // epochdays
    julian_date_at_epoch: c_double,          // jdsatepoch
    mean_motion_second_derivative: c_double, // nddot
    mean_motion_first_derivative: c_double,  // ndot

    // SGP4-type drag coefficient
    bstar: c_double,

    revolution_count_at_epoch: c_double, // rcse
    // Inclination
    pub(crate) inclination: c_double, // inclo
    // Right ascension of the ascending node
    pub(crate) raan: c_double, // nodeo
    // Eccentricity
    pub(crate) eccentricity: c_double, // ecco
    // Argument of perigee
    pub(crate) argument_of_perigee: c_double, // argpo
    // Mean anomaly
    pub(crate) mean_anomaly: c_double, // mo

    pub(crate) mean_motion: c_double, // no
}

/// Determine if two C doubles are "close", as defined by the EPSILON constant.
pub(crate) fn close(a: c_double, b: c_double) -> bool {
    let err = (a - b).abs();
    err <= EPSILON
}

impl PartialEq for OrbitalElementSet {
    fn eq(&self, other: &Self) -> bool {
        self.catalog_number == other.catalog_number
            && self.epoch_year == other.epoch_year
            && self.epochtynumrev == other.epochtynumrev
            && self.error == other.error
            && self.operation_mode == other.operation_mode
            && self.init == other.init
            && self.method == other.method
            && self.isimp == other.isimp
            && close(self.aycof, other.aycof)
            && close(self.con41, other.con41)
            && close(self.cc1, other.cc1)
            && close(self.cc4, other.cc4)
            && close(self.cc5, other.cc5)
            && close(self.d2, other.d2)
            && close(self.d3, other.d3)
            && close(self.d4, other.d4)
            && close(self.delmo, other.delmo)
            && close(self.eta, other.eta)
            && close(self.argpdot, other.argpdot)
            && close(self.omgcof, other.omgcof)
            && close(self.sinmao, other.sinmao)
            && close(self.t, other.t)
            && close(self.t2cof, other.t2cof)
            && close(self.t3cof, other.t3cof)
            && close(self.t4cof, other.t4cof)
            && close(self.t5cof, other.t5cof)
            && close(self.x1mth2, other.x1mth2)
            && close(self.x7thm1, other.x7thm1)
            && close(self.mdot, other.mdot)
            && close(self.nodedot, other.nodedot)
            && close(self.xlcof, other.xlcof)
            && close(self.xmcof, other.xmcof)
            && close(self.nodecf, other.nodecf)
            && self.irez == other.irez
            && close(self.d2201, other.d2201)
            && close(self.d2211, other.d2211)
            && close(self.d3210, other.d3210)
            && close(self.d3222, other.d3222)
            && close(self.d4410, other.d4410)
            && close(self.d4422, other.d4422)
            && close(self.d5220, other.d5220)
            && close(self.d5232, other.d5232)
            && close(self.d5421, other.d5421)
            && close(self.d5433, other.d5433)
            && close(self.dedt, other.dedt)
            && close(self.del1, other.del1)
            && close(self.del2, other.del2)
            && close(self.del3, other.del3)
            && close(self.didt, other.didt)
            && close(self.dmdt, other.dmdt)
            && close(self.dnodt, other.dnodt)
            && close(self.domdt, other.domdt)
            && close(self.e3, other.e3)
            && close(self.ee2, other.ee2)
            && close(self.peo, other.peo)
            && close(self.pgho, other.pgho)
            && close(self.pho, other.pho)
            && close(self.pinco, other.pinco)
            && close(self.plo, other.plo)
            && close(self.se2, other.se2)
            && close(self.se3, other.se3)
            && close(self.sgh2, other.sgh2)
            && close(self.sgh3, other.sgh3)
            && close(self.sgh4, other.sgh4)
            && close(self.sh2, other.sh2)
            && close(self.sh3, other.sh3)
            && close(self.si2, other.si2)
            && close(self.si3, other.si3)
            && close(self.sl2, other.sl2)
            && close(self.sl3, other.sl3)
            && close(self.sl4, other.sl4)
            && close(self.gsto, other.gsto)
            && close(self.xfact, other.xfact)
            && close(self.xgh2, other.xgh2)
            && close(self.xgh3, other.xgh3)
            && close(self.xgh4, other.xgh4)
            && close(self.xh2, other.xh2)
            && close(self.xh3, other.xh3)
            && close(self.xi2, other.xi2)
            && close(self.xi3, other.xi3)
            && close(self.xl2, other.xl2)
            && close(self.xl3, other.xl3)
            && close(self.xl4, other.xl4)
            && close(self.xlamo, other.xlamo)
            && close(self.zmol, other.zmol)
            && close(self.zmos, other.zmos)
            && close(self.atime, other.atime)
            && close(self.xli, other.xli)
            && close(self.xni, other.xni)
            && close(self.semi_major_axis, other.semi_major_axis)
            && close(self.altitude_of_periapsis, other.altitude_of_periapsis)
            && close(self.altitude_of_apoapsis, other.altitude_of_apoapsis)
            && close(self.epoch_days, other.epoch_days)
            && close(self.julian_date_at_epoch, other.julian_date_at_epoch)
            && close(
                self.mean_motion_second_derivative,
                other.mean_motion_second_derivative,
            )
            && close(
                self.mean_motion_first_derivative,
                other.mean_motion_first_derivative,
            )
            && close(self.bstar, other.bstar)
            && close(
                self.revolution_count_at_epoch,
                other.revolution_count_at_epoch,
            )
            && close(self.inclination, other.inclination)
            && close(self.raan, other.raan)
            && close(self.eccentricity, other.eccentricity)
            && close(self.argument_of_perigee, other.argument_of_perigee)
            && close(self.mean_anomaly, other.mean_anomaly)
            && close(self.mean_motion, other.mean_motion)
    }
}

impl OrbitalElementSet {
    pub fn epoch(&self) -> DateTime<Utc> {
        julian_day_to_datetime(self.julian_date_at_epoch)
    }

    pub(crate) fn mean_motion(&self) -> f64 {
        self.mean_motion as _
    }

    pub(crate) fn into_validated_result(self) -> Result<OrbitalElementSet, Error> {
        match self.error {
            0 => Ok(self),
            1 => Err(Error::InvalidEccentricity),
            2 => Err(Error::NegativeMeanMotion),
            3 => Err(Error::PertEccentricity),
            4 => Err(Error::NegativeSemiLatus),
            5 => Err(Error::SubOrbital),
            6 => Err(Error::SatelliteDecay),
            _ => Err(Error::Unknown),
        }
    }
}

pub(crate) fn julian_day_to_datetime(jd: c_double) -> DateTime<Utc> {
    let mut year = c_int::default();
    let mut month = c_int::default();
    let mut day = c_int::default();
    let mut hour = c_int::default();
    let mut minute = c_int::default();
    let mut second = c_double::default();

    unsafe {
        invjday(
            jd,
            &mut year,
            &mut month,
            &mut day,
            &mut hour,
            &mut minute,
            &mut second,
        );
    }

    let usec = (second.fract() * 1e6).floor();
    NaiveDate::from_ymd_opt(year, month as u32, day as u32)
        .unwrap()
        .and_hms_micro_opt(hour as u32, minute as u32, second as u32, usec as u32)
        .unwrap()
        .and_utc()
}

pub(crate) fn datetime_to_julian_day(d: DateTime<Utc>) -> c_double {
    let mut jd = c_double::default();

    unsafe {
        jday(
            d.year() as c_int,
            d.month() as c_int,
            d.day() as c_int,
            d.hour() as c_int,
            d.minute() as c_int,
            d.second() as c_double,
            &mut jd,
        );
    }

    jd
}

pub(crate) fn to_orbital_elements(
    line1: &str,
    line2: &str,
    rt: RunType,
    om: OperationMode,
    gc: GravitationalConstant,
) -> Result<OrbitalElementSet, Error> {
    let l1 = CString::new(line1)?;
    let l2 = CString::new(line2)?;

    let mut startmfe: c_double = 0.;
    let mut stopmfe: c_double = 0.;
    let mut deltamin: c_double = 0.;
    let mut satrec: OrbitalElementSet = OrbitalElementSet::default();

    unsafe {
        twoline2rv(
            l1.as_ptr(),
            l2.as_ptr(),
            rt.to_char(),
            rt.to_input_char(),
            om.to_char(),
            gc,
            &mut startmfe,
            &mut stopmfe,
            &mut deltamin,
            &mut satrec,
        )
    };

    satrec.into_validated_result()
}

type Vec3 = [c_double; 3];
type VectorPair = (Vec3, Vec3);

pub(crate) fn run_sgp4(
    satrec: OrbitalElementSet,
    gc: GravitationalConstant,
    min_since_epoch: f64,
) -> Result<VectorPair, Error> {
    let mut satrec_copy = satrec.to_owned();

    let mut ro: Vec3 = [0.; 3];
    let mut vo: Vec3 = [0.; 3];

    let success = unsafe {
        sgp4(
            gc,
            &mut satrec_copy,
            min_since_epoch as c_double,
            ro.as_mut_ptr(),
            vo.as_mut_ptr(),
        )
    };

    if !success {
        match satrec_copy.into_validated_result() {
            Ok(_) => unreachable!("SGP4 failure but didn't return an error"),
            Err(e) => Err(e),
        }
    } else {
        Ok((ro, vo))
    }
}

#[derive(Debug)]
pub(crate) struct ClassicalOrbitalElements {
    pub p: c_double,       // semilatus rectum               km
    pub a: c_double,       // semimajor axis                 km
    pub ecc: c_double,     // eccentricity
    pub incl: c_double,    // inclination                    0.0  to pi rad
    pub omega: c_double,   // longitude of ascending node    0.0  to 2pi rad
    pub argp: c_double,    // argument of perigee            0.0  to 2pi rad
    pub nu: c_double,      // true anomaly                   0.0  to 2pi rad
    pub m: c_double,       // mean anomaly                   0.0  to 2pi rad
    pub arglat: c_double,  // argument of latitude      (ci) 0.0  to 2pi rad
    pub truelon: c_double, // true longitude            (ce) 0.0  to 2pi rad
    pub lonper: c_double,  // longitude of periapsis    (ee) 0.0  to 2pi rad
}

#[allow(clippy::many_single_char_names)]
pub(crate) fn to_classical_elements(r: &Vec3, v: &Vec3) -> ClassicalOrbitalElements {
    let grav_consts = gravitational_constants();

    let mut p: c_double = 0.0;
    let mut a: c_double = 0.0;
    let mut ecc: c_double = 0.0;
    let mut incl: c_double = 0.0;
    let mut omega: c_double = 0.0;
    let mut argp: c_double = 0.0;
    let mut nu: c_double = 0.0;
    let mut m: c_double = 0.0;
    let mut arglat: c_double = 0.0;
    let mut truelon: c_double = 0.0;
    let mut lonper: c_double = 0.0;

    unsafe {
        rv2coe(
            r.as_ptr(),
            v.as_ptr(),
            grav_consts.mu,
            &mut p,
            &mut a,
            &mut ecc,
            &mut incl,
            &mut omega,
            &mut argp,
            &mut nu,
            &mut m,
            &mut arglat,
            &mut truelon,
            &mut lonper,
        );
    }

    ClassicalOrbitalElements {
        p,
        a,
        ecc,
        incl,
        omega,
        argp,
        nu,
        m,
        arglat,
        truelon,
        lonper,
    }
}

pub(crate) fn datetime_to_gstime(d: DateTime<Utc>) -> c_double {
    let jd = datetime_to_julian_day(d);
    unsafe { gstime(jd) }
}

pub(crate) struct GravitationalConstants {
    pub tumin: c_double,
    pub mu: c_double,
    pub radiusearthkm: c_double,
    pub xke: c_double,
    pub j2: c_double,
    pub j3: c_double,
    pub j4: c_double,
    pub j3oj2: c_double,
}

pub(crate) fn gravitational_constants() -> GravitationalConstants {
    let mut consts = GravitationalConstants {
        tumin: 0.0,
        mu: 0.0,
        radiusearthkm: 0.0,
        xke: 0.0,
        j2: 0.0,
        j3: 0.0,
        j4: 0.0,
        j3oj2: 0.0,
    };
    unsafe {
        getgravconst(
            GravitationalConstant::Wgs84,
            &mut consts.tumin,
            &mut consts.mu,
            &mut consts.radiusearthkm,
            &mut consts.xke,
            &mut consts.j2,
            &mut consts.j3,
            &mut consts.j4,
            &mut consts.j3oj2,
        );
    }

    consts
}

#[link(name = "sgp4", kind = "static")]
#[allow(non_snake_case)]
#[allow(dead_code)]
extern "C" {
    // Defined in sgp4unit.cpp

    fn sgp4init(
        whichconst: GravitationalConstant,
        opsmode: c_char,
        satn: c_int,
        epoch: c_double,
        xbstar: c_double,
        xecco: c_double,
        xargpo: c_double,
        xinclo: c_double,
        xmo: c_double,
        xno: c_double,
        xnodeo: c_double,
        satrec: &mut OrbitalElementSet,
    ) -> bool;

    fn sgp4(
        whichconst: GravitationalConstant,
        satrec: &mut OrbitalElementSet,
        tsince: c_double,
        r: *mut c_double,
        v: *mut c_double,
    ) -> bool;

    fn gstime(jdut1: c_double) -> c_double;

    fn getgravconst(
        whichconst: GravitationalConstant,
        tumin: &mut c_double,
        mu: &mut c_double,
        radiusearthkm: &mut c_double,
        xke: &mut c_double,
        j2: &mut c_double,
        j3: &mut c_double,
        j4: &mut c_double,
        j3oj2: &mut c_double,
    );

    // Defined in sgp4io.cpp

    fn twoline2rv(
        longstr1: *const c_char,
        longstr2: *const c_char,
        typerun: c_char,
        typeinput: c_char,
        opsmode: c_char,
        whichconst: GravitationalConstant,
        startmfe: &mut c_double,
        stopmfe: &mut c_double,
        deltamin: &mut c_double,
        satrec: &mut OrbitalElementSet,
    );

    // Defined in sgp4ext.cpp

    fn sgn(x: c_double) -> c_double;

    fn mag(x: *const c_double) -> c_double;

    fn cross(v1: *const c_double, v2: *const c_double, out: *mut c_double);

    fn dot(v1: *const c_double, v2: *const c_double) -> c_double;

    fn angle(v1: *const c_double, v2: *const c_double) -> c_double;

    fn newtonnu(ecc: c_double, nu: c_double, e0: &mut c_double, m: &mut c_double);

    fn asinh(xval: c_double) -> c_double;

    fn rv2coe(
        r: *const c_double, // [c_double; 3]
        v: *const c_double, // [c_double; 3]
        mu: c_double,
        p: &mut c_double,
        a: &mut c_double,
        ecc: &mut c_double,
        incl: &mut c_double,
        omega: &mut c_double,
        argp: &mut c_double,
        nu: &mut c_double,
        m: &mut c_double,
        arglat: &mut c_double,
        truelon: &mut c_double,
        lonper: &mut c_double,
    );

    fn jday(
        year: c_int,
        mon: c_int,
        day: c_int,
        hr: c_int,
        minute: c_int,
        sec: c_double,
        jd: &mut c_double,
    );

    fn days2mdhms(
        year: c_int,
        days: c_double,
        mon: &mut c_int,
        day: &mut c_int,
        hr: &mut c_int,
        minute: &mut c_int,
        sec: &mut c_double,
    );

    fn invjday(
        julian_day: c_double,
        year: &mut c_int,
        month: &mut c_int,
        day: &mut c_int,
        hour: &mut c_int,
        minute: &mut c_int,
        seconds: &mut c_double,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_propagation() {
        // twoline2rv( longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
        //             startmfe, stopmfe, deltamin, satrec );
        // fprintf(outfile, "%ld xx\n", satrec.satnum);
        // printf(" %ld\n", satrec.satnum);
        // // call the propagator to get the initial state vector value
        // sgp4 (whichconst, satrec,  0.0, ro,  vo);

        // ISS (ZARYA)
        // 1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992
        // 2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767

        let longstr1 =
            CString::new("1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992")
                .unwrap();
        let longstr2 =
            CString::new("2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767")
                .unwrap();
        let typerun = 'v' as c_char;
        let typeinput = 'v' as c_char;
        let opsmode = 'a' as c_char;
        let whichconst = GravitationalConstant::Wgs72;
        let mut startmfe: c_double = 0.;
        let mut stopmfe: c_double = 0.;
        let mut deltamin: c_double = 0.;
        let mut satrec: OrbitalElementSet = OrbitalElementSet::default();
        let mut ro: [c_double; 3] = [0.; 3];
        let mut vo: [c_double; 3] = [0.; 3];

        unsafe {
            twoline2rv(
                longstr1.as_ptr(),
                longstr2.as_ptr(),
                typerun,
                typeinput,
                opsmode,
                whichconst,
                &mut startmfe,
                &mut stopmfe,
                &mut deltamin,
                &mut satrec,
            );
        }

        let satrec = satrec;
        let mut satrec_copy = satrec.to_owned();

        unsafe {
            sgp4(
                whichconst,
                &mut satrec_copy,
                0.0,
                ro.as_mut_ptr(),
                vo.as_mut_ptr(),
            );
        }

        assert_eq!(satrec_copy.error, 0);

        // Propagate out 60 minutes.
        let mut satrec_copy = satrec.to_owned();

        unsafe {
            sgp4(
                whichconst,
                &mut satrec_copy,
                60.0,
                ro.as_mut_ptr(),
                vo.as_mut_ptr(),
            );
        }
        assert_eq!(satrec_copy.error, 0);
    }

    #[test]
    fn test_errors() {
        //This TLE decays on by the end Feb 2022
        // 0 LEMUR 2 MCCULLAGH
        //1 43051U 17071Q   22046.92182028  .07161566  12340-4  74927-3 0  9993
        //2 43051  51.6207 236.5853 0009084 284.2762  75.7254 16.36736354237455
        let longstr1 =
            CString::new("1 43051U 17071Q   22046.92182028  .07161566  12340-4  74927-3 0  9993")
                .unwrap();
        let longstr2 =
            CString::new("2 43051  51.6207 236.5853 0009084 284.2762  75.7254 16.36736354237455")
                .unwrap();
        let typerun = 'v' as c_char;
        let typeinput = 'v' as c_char;
        let opsmode = 'a' as c_char;
        let whichconst = GravitationalConstant::Wgs72;
        let mut startmfe: c_double = 0.;
        let mut stopmfe: c_double = 0.;
        let mut deltamin: c_double = 0.;
        let mut satrec: OrbitalElementSet = OrbitalElementSet::default();
        let mut ro: [c_double; 3] = [0.; 3];
        let mut vo: [c_double; 3] = [0.; 3];

        unsafe {
            twoline2rv(
                longstr1.as_ptr(),
                longstr2.as_ptr(),
                typerun,
                typeinput,
                opsmode,
                whichconst,
                &mut startmfe,
                &mut stopmfe,
                &mut deltamin,
                &mut satrec,
            );
        }

        let satrec = satrec;
        let mut satrec_copy = satrec.to_owned();

        unsafe {
            sgp4(
                whichconst,
                &mut satrec_copy,
                0.0,
                ro.as_mut_ptr(),
                vo.as_mut_ptr(),
            );
        }

        assert_eq!(satrec_copy.error, 0);

        // Propagate out 60 minutes.
        let mut satrec_copy = satrec.to_owned();
        unsafe {
            sgp4(
                whichconst,
                &mut satrec_copy,
                10000.,
                ro.as_mut_ptr(),
                vo.as_mut_ptr(),
            );
        }
        assert_eq!(satrec_copy.error, 6)
    }

    #[test]
    fn test_close() {
        assert!(!close(0.1, 0.3));
        assert!(close(0.5, 0.5));
        // Note: Due to floating point errors, these tests can't use 0.5 - EPSILON. We exploit
        // argument ordering here to make sure order doesn't matter while staying on the boundary of
        // the acceptable region.
        assert!(close(0.5, 0.5 + EPSILON));
        assert!(close(0.5 + EPSILON, 0.5));
    }
}
