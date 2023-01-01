const std = @import("std");
const Builder = std.build.Builder;
const Pkg = std.build.Pkg;

const pkgs = struct {
    const zstats = Pkg{
        .name = "zstats",
        .source = .{ .path = "./zstats.zig" },
        .dependencies = &[_]Pkg{ },
    };
};

pub fn build(b: *std.build.Builder) void {
    const target = b.standardTargetOptions(.{});
    const mode = b.standardReleaseOptions();

    const exe_targets = [_]Target{
        // .{ .name = "foo", .src = "foo.zig", .desc ="foo"},
    };

    for (exe_targets) |e_target| {
        e_target.build(b);
    }

    const exe_tests = b.addTest("./zstats.zig");
    exe_tests.setTarget(target);
    exe_tests.setBuildMode(mode);

    exe_tests.addPackage(pkgs.zstats);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&exe_tests.step);
}

//---------------------------------------------------------------------------

const Target = struct {
    name: []const u8,
    src: []const u8,
    desc: []const u8,

    pub fn build(self: Target, b: *std.build.Builder) void {

        // write executables here instead
        b.exe_dir = b.pathFromRoot("./bin");

        // const target = b.standardTargetOptions(.{});
        const mode = b.standardReleaseOptions();

        var exe = b.addExecutable(self.name, self.src);
        exe.setBuildMode(b.standardReleaseOptions());
        // exe.setTarget(target);
        exe.setBuildMode(mode);

        exe.addPackage(pkgs.zstats);

        exe.install();

        const run_cmd = exe.run();
        run_cmd.step.dependOn(b.getInstallStep());
        if (b.args) |args| {
            run_cmd.addArgs(args);
        }
    }
};
