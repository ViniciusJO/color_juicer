const std = @import("std");

test "COLORS" {
    var cor = RGBA.init(.{ .hex = 0x12345678 });
    std.debug.print("\n\nHEX: #{X}\nRGB: {}\nLCH: {}\nRGB: {}(converted back)\n\n", .{ cor.to_u32(), cor, cor.to_lch(), cor.to_lch().to_rgba() });

    cor = RGBA.init(.{ .hex = 0xFF0000FF });
    std.debug.print("\n\nHEX: #{X}\nRGB: {}\nLCH: {}\nRGB: {}(converted back)\n\n", .{ cor.to_u32(), cor, cor.to_lch(), cor.to_lch().to_rgba() });
}

pub const RGBA = struct {
    r: u8,
    g: u8,
    b: u8,
    a: u8,

    const Self = @This();

    pub fn init(color: AnyColor) Self {
        return switch(color) {
            .rgba => |c| c,
            .lch  => |c| RGBA.from_lch(c),
            .hex  => |c| Self{
                .r = @intCast((c >> 24)&0xFF),
                .g = @intCast((c >> 16)&0xFF),
                .b = @intCast((c >>  8)&0xFF),
                .a = @intCast((c >>  0)&0xFF),
            },
        };
    }

    pub fn to_lch(self: *Self) LCH {
        const rf = u8_to_f32(self.r);
        const gf = u8_to_f32(self.g);
        const bf = u8_to_f32(self.b);

        const rl = srgb_to_linear(rf);
        const gl = srgb_to_linear(gf);
        const bl = srgb_to_linear(bf);

        const l = 0.4122214708 * rl + 0.5363325363 * gl + 0.0514459929 * bl;
        const m = 0.2119034982 * rl + 0.6806995451 * gl + 0.1073969566 * bl;
        const s = 0.0883024619 * rl + 0.2817188376 * gl + 0.6299787005 * bl;

        const l_ = std.math.cbrt(l);
        const m_ = std.math.cbrt(m);
        const s_ = std.math.cbrt(s);

        const L = 0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_;
        const A = 1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_;
        const B = 0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_;

        return .{
            .l = L,
            .c = @sqrt(A * A + B * B),
            .h = std.math.atan2(B, A)/std.math.pi*180,
            .a = u8_to_f32(self.a),
        };
    }

    pub fn from_lch(lch: LCH) Self {
        const h = lch.h/180*std.math.pi;
        const A = lch.c * @cos(h);
        const B = lch.c * @sin(h);

        const l_ = lch.l + 0.3963377774 * A + 0.2158037573 * B;
        const m_ = lch.l - 0.1055613458 * A - 0.0638541728 * B;
        const s_ = lch.l - 0.0894841775 * A - 1.2914855480 * B;

        const l = l_ * l_ * l_;
        const m = m_ * m_ * m_;
        const s = s_ * s_ * s_;

        const r_lin =  4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
        const g_lin = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
        const b_lin = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s;

        return .{
            .r = f32_to_u8(linear_to_srgb(r_lin)),
            .g = f32_to_u8(linear_to_srgb(g_lin)),
            .b = f32_to_u8(linear_to_srgb(b_lin)),
            .a = f32_to_u8(lch.a),
        };
    }

    pub fn to_hex(self: *Self) []const u8 {
        var str = [1]u8{0}**9;
        _ = std.fmt.bufPrint(&str, "#{:0>2}{:0>2}{:0>2}{:0>2}", .{ self.r, self.g, self.b, self.a }) catch return "";
        return &str;
    }

    pub fn to_u32(self: *Self) u32 {
        return
            @as(u32, @intCast(self.r)) << 24 |
            @as(u32, @intCast(self.g)) << 16 |
            @as(u32, @intCast(self.b)) <<  8 |
            @as(u32, @intCast(self.a)) <<  0;
    }
};

pub const LCH = struct {
    l: f32,
    c: f32,
    h: f32, // radians
    a: f32, // alpha 0..1
    
    const Self = @This();

    pub fn init(color: AnyColor) Self {
        return switch(color) {
            .rgba => |c| c.to_lch(),
            .lch  => |c| c,
            .hex  => |c| RGBA.init(c).to_lch(),
        };
    }

    pub fn to_rgba(self: *const Self) RGBA {
        return RGBA.from_lch(self.*);
    }

    pub fn from_rgba(rgba: RGBA) Self {
        return rgba.to_lch();
    }
};

pub const AnyColor = union(enum){ rgba: RGBA, lch: LCH, hex: u32 };

// ==================== UTILS ====================

fn clamp01(x: f32) f32 {
    return if (x < 0) 0 else if (x > 1) 1 else x;
}

fn u8_to_f32(x: u8) f32 {
    return @as(f32, @floatFromInt(x)) / 255.0;
}

fn f32_to_u8(x: f32) u8 {
    return @intFromFloat(@round(clamp01(x) * 255.0));
}

fn srgb_to_linear(c: f32) f32 {
    return if (c <= 0.04045)
        c / 12.92
    else
        std.math.pow(f32, (c + 0.055) / 1.055, 2.4);
}

fn linear_to_srgb(c: f32) f32 {
    return if (c <= 0.0031308)
        12.92 * c
    else
        1.055 * std.math.pow(f32, c, 1.0 / 2.4) - 0.055;
}


