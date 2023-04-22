extern crate cc;

fn main() {
    if cfg!(feature = "latest-gcc") {
        if cfg!(feature = "blas") {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/part_linpack.f")
                .file("src/opkda2.f")
                .file("src/opkda1.f")
                .file("src/opkdmain.f")
                .compile("libodepack.a");
            println!("cargo:rustc-link-lib=dylib=blas");
        } else {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/part_blas.f")
                .file("src/part_linpack.f")
                .file("src/opkda2.f")
                .file("src/opkda1.f")
                .file("src/opkdmain.f")
                .compile("libodepack.a");
        }
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                // use for version of gfortran 10+
                .flag("-fallow-argument-mismatch")
                .file("src/dc_lapack.f")
                .file("src/radau.f")
                .compile("libradau.a");
            println!("cargo:rustc-link-lib=dylib=lapack");
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
        if cfg!(feature = "blas") {
            cc::Build::new()
                .flag("-w")
                .file("src/part_linpack.f")
                .file("src/opkda2.f")
                .file("src/opkda1.f")
                .file("src/opkdmain.f")
                .compile("libodepack.a");
        } else {
            cc::Build::new()
                .flag("-w")
                .file("src/part_blas.f")
                .file("src/part_linpack.f")
                .file("src/opkda2.f")
                .file("src/opkda1.f")
                .file("src/opkdmain.f")
                .compile("libodepack.a");
        }
        if cfg!(feature = "lapack") {
            cc::Build::new()
                .flag("-w")
                .file("src/dc_lapack.f")
                .file("src/radau.f")
                .compile("libradau.a");
            println!("cargo:rustc-link-lib=dylib=lapack");
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
