extern crate cc;

fn main() {
    if cfg!(feature = "latest-gcc") {
        cc::Build::new()
            .flag("-w")
            // use for version of gfortran 10+
            .flag("-fallow-argument-mismatch")
            .file("src/odepack.f")
            .compile("libodepack.a");
    } else {
        cc::Build::new()
            .flag("-w")
            .file("src/odepack.f")
            .compile("libodepack.a");
    }

    println!("cargo:rustc-link-lib=static=odepack");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
