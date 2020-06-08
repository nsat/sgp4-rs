//! # SGP-4 System
//!
//! This module contains FFI bindings to the Vallado C++ implementation of SGP-4. These bindings
//! rely heavily on `unsafe` code, and it is recommended to instead use the other, high level
//! bindings provided by the `sgp4` crate.

extern crate libc;

use libc::{c_char, c_double, c_int, c_long};
use std::ffi::{CString, NulError};

use chrono::prelude::*;
use chrono::DateTime;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    CStringNulError(#[from] NulError),
    #[error("Failed to convert two-line element to orbital element set")]
    TwoLine2Rv,
    #[error("Failure in SGP4 propagator")]
    Sgp4,
}

#[allow(dead_code)]
pub enum RunType {
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

impl Default for RunType {
    fn default() -> Self {
        RunType::Verification
    }
}

#[allow(dead_code)]
pub enum InputType {
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
pub enum OperationMode {
    AirForceSpaceCenter,
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

impl Default for OperationMode {
    fn default() -> Self {
        OperationMode::Improved
    }
}

const EPSILON: f64 = 0.000_000_1;

#[repr(C)]
#[allow(non_camel_case_types)]
#[allow(dead_code)]
#[derive(Copy, Clone, Debug)]
pub enum GravitationalConstant {
    Wgs72Old,
    Wgs72,
    Wgs84,
}

#[repr(C)]
#[allow(dead_code)]
#[derive(Default, Clone, Copy, Debug)]
pub struct OrbitalElementSet {
    satnum: c_long,
    epochyr: c_int,
    epochtynumrev: c_int,

    error: c_int,

    operationmode: c_char,
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

    a: c_double,
    altp: c_double,
    alta: c_double,
    epochdays: c_double,
    jdsatepoch: c_double,
    nddot: c_double,
    ndot: c_double,

    bstar: c_double,
    rcse: c_double,
    inclo: c_double,
    nodeo: c_double,
    ecco: c_double,
    argpo: c_double,
    mo: c_double,

    no: c_double,
}

/// Determine if two C doubles are "close", as defined by the EPSILON constant.
fn close(a: c_double, b: c_double) -> bool {
    let err = (a - b as f64).abs();
    err <= EPSILON
}

impl PartialEq for OrbitalElementSet {
    fn eq(&self, other: &Self) -> bool {
        self.satnum == other.satnum
            && self.epochyr == other.epochyr
            && self.epochtynumrev == other.epochtynumrev
            && self.error == other.error
            && self.operationmode == other.operationmode
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
            && close(self.a, other.a)
            && close(self.altp, other.altp)
            && close(self.alta, other.alta)
            && close(self.epochdays, other.epochdays)
            && close(self.jdsatepoch, other.jdsatepoch)
            && close(self.nddot, other.nddot)
            && close(self.ndot, other.ndot)
            && close(self.bstar, other.bstar)
            && close(self.rcse, other.rcse)
            && close(self.inclo, other.inclo)
            && close(self.nodeo, other.nodeo)
            && close(self.ecco, other.ecco)
            && close(self.argpo, other.argpo)
            && close(self.mo, other.mo)
            && close(self.no, other.no)
    }
}

impl OrbitalElementSet {
    pub fn epoch(&self) -> DateTime<Utc> {
        let mut year = c_int::default();
        let mut month = c_int::default();
        let mut day = c_int::default();
        let mut hour = c_int::default();
        let mut minute = c_int::default();
        let mut second = c_double::default();

        unsafe {
            invjday(
                self.jdsatepoch,
                &mut year,
                &mut month,
                &mut day,
                &mut hour,
                &mut minute,
                &mut second,
            );
        }

        Utc.ymd(year, month as u32, day as u32).and_hms(
            hour as u32,
            minute as u32,
            second as u32,
        )
    }
}

pub fn to_orbital_elements(
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

    match satrec.error {
        0 => Ok(satrec),
        // TODO Expand this match to include specific error conditions
        _ => Err(Error::TwoLine2Rv),
    }
}

type VectorPair = ([c_double; 3], [c_double; 3]);

pub fn run_sgp4(
    satrec: OrbitalElementSet,
    gc: GravitationalConstant,
    min_since_epoch: f64,
) -> Result<VectorPair, Error> {
    let mut satrec_copy = satrec.to_owned();

    let mut ro: [c_double; 3] = [0.; 3];
    let mut vo: [c_double; 3] = [0.; 3];

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
        Err(Error::Sgp4)
    } else {
        Ok((ro, vo))
    }
}

#[link(name = "sgp4", kind = "static")]
#[allow(non_snake_case)]
#[allow(dead_code)]
extern "C" {
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
