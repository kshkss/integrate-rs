extern crate cc;

fn main() {
    if cfg!(feature = "latest-gcc") {
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/odepack_sub2.f")
            .file("src/odepack_sub1.f")
            .file("src/odepack.f")
            .compile("libodepack.a");
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/dc_lapack.f")
                .file("src/radau.f")
                .compile("libradau.a");
        } else {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/decsol.f")
                .file("src/dc_decsol.f")
                .file("src/radau.f")
                .compile("libradau.a");
        }
    } else {
        cc::Build::new()
            .flag("-w")
            .file("src/odepack_sub2.f")
            .file("src/odepack_sub1.f")
            .file("src/odepack.f")
            .compile("libodepack.a");
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                .file("src/dc_lapack.f")
                .file("src/radau.f")
                .compile("libradau.a");
        } else {
            cc::Build::new()
                .flag("-w")
                .file("src/decsol.f")
                .file("src/dc_decsol.f")
                .file("src/radau.f")
                .compile("libradau.a");
        }
    }

    println!("cargo:rustc-link-lib=static=odepack");
    println!("cargo:rustc-link-lib=static=radau");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
