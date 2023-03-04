extern crate cc;

fn main() {
    if cfg!(feature = "latest-gcc") {
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/odepack.f")
            .compile("libodepack.a");
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/radau.f")
            .compile("libradau.a");
    } else {
        cc::Build::new()
            .flag("-w")
            .file("src/odepack.f")
            .compile("libodepack.a");
        cc::Build::new()
            .flag("-w")
            .file("src/radau.f")
            .compile("libradau.a");
    }

    println!("cargo:rustc-link-lib=static=odepack");
    println!("cargo:rustc-link-lib=static=radau");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
