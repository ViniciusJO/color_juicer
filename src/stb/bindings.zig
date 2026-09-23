const Self = @This();

pub const image = struct {
    pub const load = @extern(
        *const fn(filename: [*c]const u8, x: [*c]c_int, y: [*c]c_int, channels_in_file: [*c]c_int, desired_channels: c_int) callconv(.c) [*c]u8,
        .{ .name = "stbi_load" }
    );

    const ImageError = error {
        LoadFailed,
        WriteFailed,
    };

    const Image = struct {
        path: []const u8,
        width: u32,
        height: u32,
        channels: u8,
        pixels: []const u8,
    };

    const Format = enum { png, bmp, tga, hdr, jpg, };

    pub fn load_(filename: []const u8, desired_channels: u8) ImageError!Image {
        var x: c_int = 0;
        var y: c_int = 0;
        var ch: c_int = 0;
        const res = load(filename.ptr, &x, &y, &ch, @intCast(desired_channels));
        if(res == null) return ImageError.LoadFailed;
        return Image{
            .path = filename,
            .width = @intCast(x),
            .height= @intCast(y),
            .channels = @intCast(ch),
            .pixels = res[0..@as(usize, @intCast(x*y*ch))],
        };
    }

    pub const write_png = @extern(
        *const fn(filename: [*c]const u8, w: c_int, h: c_int, comp: c_int, data: ?*const anyopaque, stride_in_bytes: c_int) callconv(.c) c_int,
        .{ .name = "stbi_write_png" },
    );

    pub fn write(img: Image, filename: []const u8, format: Format) ImageError!void {
        const res = switch(format) {
            .png => write_png(filename, img.width, img.height, img.channels, img.pixels.ptr, img.width*img.channels*@sizeOf(img.pixels[0])),
            else => 0,
        };
        if(res == 0) return ImageError.WriteFailed;
    }

    pub const resize_uint8_linear = @extern(
        *const fn(input_pixels: [*c]const u8, input_w: c_int, input_h: c_int, input_stride_in_bytes: c_int, output_pixels: [*c]u8, output_w: c_int, output_h: c_int, output_stride_in_bytes: c_int, pixel_type: c_uint) callconv(.c) [*c]u8,
        .{ .name = "stbir_resize_uint8_linear" },
    );
};
