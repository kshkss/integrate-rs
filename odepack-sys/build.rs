extern crate cc;

fn main() {
    if cfg!(feature = "blas") {
        cc::Build::new()
            .flag("-w")
            .file("src/blas.f")
            .compile("libblas.a");
        println!("cargo:rustc-link-lib=static=blas");
    } else {
        println!("cargo:rustc-link-lib=dylib=blas");
    }
    if cfg!(feature = "linpack") {
        cc::Build::new()
            .flag("-w")
            .file("src/linpack.f")
            .compile("liblinpack.a");
        println!("cargo:rustc-link-lib=static=linpack");
    } else {
        println!("cargo:rustc-link-lib=dylib=linpack");
    }
    if cfg!(feature = "slatec") {
        cc::Build::new()
            .flag("-w")
            .file("src/slatec.f")
            .compile("libslatec.a");
        println!("cargo:rustc-link-lib=static=slatec");
    } else {
        println!("cargo:rustc-link-lib=dylib=slatec");
    }
    if cfg!(feature = "latest-gcc") {
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/opkda1.f")
            .file("src/opkdmain.f")
            .compile("libodepack.a");
    } else {
        cc::Build::new()
            .flag("-w")
            .file("src/opkda1.f")
            .file("src/opkdmain.f")
            .compile("libodepack.a");
    }

    println!("cargo:rustc-link-lib=static=odepack");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
