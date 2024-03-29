const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const unit_tests = b.addTest(.{
        .root_source_file = .{.path = "./zstats.zig"},
        .target = target,
        .optimize = optimize,
    });

    unit_tests.addAnonymousModule("zstats", .{
        .source_file = .{.path = "./zstats.zig"},
    });

    const unit_tests_step = b.step("test", "Run unit tests");
    unit_tests_step.dependOn(&unit_tests.step);

}

