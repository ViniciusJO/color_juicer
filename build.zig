const std = @import("std");

pub fn build(b: *std.Build) void {
    b.verbose = true;
    // b.debug_incremental = true;

    const target = b.resolveTargetQuery(.{ .cpu_arch = .x86_64, .os_tag = .linux, .cpu_model = .baseline });
    const optimize = b.standardOptimizeOption(.{});

    const wf = b.addWriteFiles();
    const stb_imp = wf.add("stb_implementations.c",
        \\#define STB_IMAGE_IMPLEMENTATION 
        \\#define STB_IMAGE_RESIZE_IMPLEMENTATION 
        \\#define STB_IMAGE_WRITE_IMPLEMENTATION 
        \\#include "../../../src/stb/stb_image.h"
        \\#include "../../../src/stb/stb_image_resize2.h"
        \\#include "../../../src/stb/stb_image_write.h"
    );

    const stb_image = b.addTranslateC(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/stb/stb_image.h"),
    });

    const stb_image_write = b.addTranslateC(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/stb/stb_image_write.h"),
    });

    const stb_image_resize = b.addTranslateC(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/stb/stb_image_resize2.h"),
    });

    const args_dep = b.dependency("args", .{ .target = target, .optimize = optimize }); //.module("args")
                                                                                        //
    const kmeans = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/main.zig"),
        .imports = &.{
            .{ .name = "args_parser", .module = args_dep.module("args") },
            .{ .name = "stbi",  .module = stb_image.createModule() },
            .{ .name = "stbiw", .module = stb_image_write.createModule() },
            .{ .name = "stbir", .module = stb_image_resize.createModule() },
        },
    });


    kmeans.addCSourceFile(.{ .file = stb_imp, .flags = &[_][]const u8{ "-fno-sanitize=undefined" } });

    // kmeans.addIncludeDir("src/stb");

    // kmeans.addCSourceFiles(.{
    //     .files = &[_][]const u8 {
    //         "src/stb/stb_image.h",
    //         "src/stb/stb_image_resize2.h",
    //         "src/stb/stb_image_write.h",
    //     },
    //     .flags = &[_][]const u8 {
    //         "-g",
    //         "-DSTB_IMAGE_IMPLEMENTATION",
    //         "-DSTB_IMAGE_WRITE_IMPLEMENTATION",
    //         "-DSTB_IMAGE_RESIZE_IMPLEMENTATION",
    //     },
    //     .language = .c
    // });

    const exe = b.addExecutable(.{
        .name = "color_juicer",
        .root_module = kmeans,
        .linkage = .static,
    });

    exe.step.dependOn(&wf.step);
    exe.step.dependOn(&stb_image.step);
    exe.step.dependOn(&stb_image_write.step);
    exe.step.dependOn(&stb_image_resize.step);

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| { run_cmd.addArgs(args); }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // const exe_unit_tests = b.addTest(.{
    //     .root_source_file = b.path("src/main.zig"),
    //     .target = target,
    //     .optimize = optimize,
    // });
    //
    // const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    // Similar to creating the run step earlier, this exposes a `test` step to
    // the `zig build --help` menu, providing a way for the user to request
    // running the unit tests.
    // const test_step = b.step("test", "Run unit tests");
    // test_step.dependOn(&run_exe_unit_tests.step);

    // const zigimg_dependency = b.dependency("zigimg", .{
    //     .target = target,
    //     .optimize = optimize,
    // });
    //
    // exe.root_module.addImport("zigimg", zigimg_dependency.module("zigimg"));
}
