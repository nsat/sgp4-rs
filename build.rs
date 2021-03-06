fn main() {
    let env_path = option_env!("SGP4_LIB_DIR");

    if let Some(sgp4_path) = env_path {
        println!("cargo:rustc-link-search=native={}", sgp4_path);
    } else {
        cc::Build::new()
            .cpp(true)
            .file("src/sgp4/sgp4ext.cpp")
            .file("src/sgp4/sgp4unit.cpp")
            .file("src/sgp4/sgp4io.cpp")
            .cpp_link_stdlib(None)
            // Uncomment to enable strict compilation.
            // This is not typically enabled because it may cause compilation failures in
            // environments we can't test, so it is primarily useful for developers of the crate.
            // .warnings(true)
            // .warnings_into_errors(true)
            .compile("libsgp4.a");
    }

    println!("cargo:rustc-link-lib=static=sgp4");
}
