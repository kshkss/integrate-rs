extern crate cc;

fn main() {
    if cfg!(feature = "latest-gcc") {
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/lapack.f")
                .file("src/lapackc.f")
                .compile("liblapack.a");
            println!("cargo:rustc-link-lib=static=lapack");
        } else {
            println!("cargo:rustc-link-lib=dylib=lapack");
        }
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/dc_lapack.f")
            .file("src/radau.f")
            .compile("libradau.a");
    } else {
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                .file("src/lapack.f")
                .file("src/lapackc.f")
                .compile("liblapack.a");
            println!("cargo:rustc-link-lib=static=lapack");
        } else {
            println!("cargo:rustc-link-lib=dylib=lapack");
        }
        cc::Build::new()
            .flag("-w")
            .file("src/dc_lapack.f")
            .file("src/radau.f")
            .compile("libradau.a");
    }

    println!("cargo:rustc-link-lib=static=radau");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
