const std = @import("std");
const math = std.math;
const Complex = std.math.Complex;

const ElementType = @import("./helpers.zig").ElementType;
const ValueType = @import("./helpers.zig").ValueType;

pub fn mean(comptime T: type, x: T) ElementType(T) {
    comptime var T_elem = ElementType(T);
    comptime var R = ValueType(T);

    var mu: T_elem = undefined;

    switch (T) {
        []f32, []f64 => {
            var sum: T_elem = 0.0;
            for (x) |xval| {
                sum += xval;
            }
            mu = sum / @as(R, @floatFromInt(x.len));
        },

        []Complex(f32), []Complex(f64) => {
            var sum = T_elem.init(0, 0);

            for (x) |xval| {
                sum.re += xval.re;
                sum.im += xval.im;
            }
            mu.re = sum.re / @as(R, @floatFromInt(x.len));
            mu.im = sum.im / @as(R, @floatFromInt(x.len));
        },
        else => {
            @compileError("type not implemented");
        },
    }
    return mu;
}

pub fn variance(comptime T: type, x: T) ValueType(T) {
    comptime var T_elem = ElementType(T);
    comptime var R = ValueType(T);

    switch (T) {
        []f32, []f64 => {
            var sum1: T_elem = 0.0;
            var sum2: T_elem = 0.0;
            for (x) |xval| {
                sum1 += xval;
                sum2 += xval * xval;
            }
            var m1: T_elem = sum1 / @as(R, @floatFromInt(x.len));
            var m2: T_elem = sum2 / @as(R, @floatFromInt(x.len));

            return m2 - m1 * m1;
        },

        []Complex(f32), []Complex(f64) => {
            var sum1 = T_elem.init(0, 0);
            var sum2: R = 0;

            var m1: T_elem = undefined;
            var m2: R = undefined;

            for (x) |xval| {
                sum1.re += xval.re;
                sum1.im += xval.im;
                sum2 += (xval.re * xval.re) + (xval.im * xval.im);
            }

            m1.re = sum1.re / @as(R, @floatFromInt(x.len));
            m1.im = sum1.im / @as(R, @floatFromInt(x.len));
            m2 = sum2 / @as(R, @floatFromInt(x.len));

            return m2 - ((m1.re * m1.re) + (m1.im * m1.im));
        },
        else => {
            @compileError("type not implemented");
        },
    }
}

pub fn mean_variance(comptime T: type, x: T, mu: *ElementType(T), sigma_sqr: *ValueType(T)) void {
    comptime var T_elem = ElementType(T);
    comptime var R = ValueType(T);

    switch (T) {
        []f32, []f64 => {
            var sum1: T_elem = 0.0;
            var sum2: T_elem = 0.0;
            for (x) |xval| {
                sum1 += xval;
                sum2 += xval * xval;
            }
            var mu_tmp: T_elem = sum1 / @as(R, @floatFromInt(x.len));
            var m2: T_elem = sum2 / @as(R, @floatFromInt(x.len));

            mu.* = mu_tmp;
            sigma_sqr.* = m2 - mu_tmp * mu_tmp;
        },

        []Complex(f32), []Complex(f64) => {
            var sum1 = T_elem.init(0, 0);
            var sum2: R = 0;

            var mu_tmp: T_elem = undefined;
            var m2: R = undefined;

            for (x) |xval| {
                sum1.re += xval.re;
                sum1.im += xval.im;
                sum2 += (xval.re * xval.re) + (xval.im * xval.im);
            }

            mu_tmp.re = sum1.re / @as(R, @floatFromInt(x.len));
            mu_tmp.im = sum1.im / @as(R, @floatFromInt(x.len));
            m2 = sum2 / @as(R, @floatFromInt(x.len));

            mu.* = mu_tmp;
            sigma_sqr.* = m2 - ((mu_tmp.re * mu_tmp.re) + (mu_tmp.im * mu_tmp.im));
        },
        else => {
            @compileError("type not implemented");
        },
    }
}

const print = std.debug.print;
const Allocator = std.mem.Allocator;
const eps = 1.0e-5;
const n = 4;

test "\t mean \t  real array\n" {
    inline for (.{ f32, f64 }) |T| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(T, n);
        x[0] = 1.0;
        x[1] = -2.0;
        x[2] = 3.0;
        x[3] = 4.0;

        var mu: T = mean([]T, x);
        try std.testing.expectApproxEqAbs(@as(T, 1.5), mu, eps);
    }
}

test "\t mean \t  complex array\n" {
    inline for (.{ f32, f64 }) |T| {
        const C: type = Complex(T);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(C, n);
        x[0].re = 1.0;
        x[1].re = -2.0;
        x[2].re = 3.0;
        x[3].re = 4.0;

        x[0].im = 2.0;
        x[1].im = 4.0;
        x[2].im = 6.0;
        x[3].im = 8.0;

        var mu: C = mean([]C, x);
        try std.testing.expectApproxEqAbs(@as(T, 1.5), mu.re, eps);
        try std.testing.expectApproxEqAbs(@as(T, 5.0), mu.im, eps);
    }
}

test "\t variance \t  real array\n" {
    inline for (.{ f32, f64 }) |T| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(T, n);
        x[0] = 1.0;
        x[1] = -2.0;
        x[2] = 3.0;
        x[3] = 4.0;

        var v: T = variance([]T, x);
        try std.testing.expectApproxEqAbs(@as(T, 5.25), v, eps);
    }
}

test "\t variance \t  complex array\n" {
    inline for (.{ f32, f64 }) |T| {
        const C = Complex(T);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(C, n);
        x[0].re = 1.0;
        x[1].re = -2.0;
        x[2].re = 3.0;
        x[3].re = 4.0;

        x[0].im = 1.0;
        x[1].im = 2.0;
        x[2].im = -2.0;
        x[3].im = 4.0;

        var v: T = variance([]C, x);
        try std.testing.expectApproxEqAbs(@as(T, 9.9375), v, eps);
    }
}

test "\t mean_variance \t real array\n" {
    inline for (.{ f32, f64 }) |T| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(T, n);
        x[0] = 1.0;
        x[1] = -2.0;
        x[2] = 3.0;
        x[3] = 4.0;

        var mu: T = undefined;
        var sigma_sqr: T = undefined;

        mean_variance([]T, x, &mu, &sigma_sqr);

        try std.testing.expectApproxEqAbs(@as(T, 1.5), mu, eps);
        try std.testing.expectApproxEqAbs(@as(T, 5.25), sigma_sqr, eps);
    }
}

test "\t mean_variance \t complex array\n" {
    inline for (.{ f32, f64 }) |T| {
        const C = Complex(T);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try allocator.alloc(C, n);
        x[0].re = 1.0;
        x[1].re = -2.0;
        x[2].re = 3.0;
        x[3].re = 4.0;

        x[0].im = 1.0;
        x[1].im = 2.0;
        x[2].im = -2.0;
        x[3].im = 4.0;

        var mu: C = undefined;
        var sigma_sqr: T = undefined;

        mean_variance([]C, x, &mu, &sigma_sqr);
        try std.testing.expectApproxEqAbs(@as(T, 1.5), mu.re, eps);
        try std.testing.expectApproxEqAbs(@as(T, 1.25), mu.im, eps);
        try std.testing.expectApproxEqAbs(@as(T, 9.9375), sigma_sqr, eps);
    }
}
