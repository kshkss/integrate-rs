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
    }

    println!("cargo:rustc-link-lib=static=odepack");
    println!("cargo:rustc-link-lib=dylib=gfortran");
}
