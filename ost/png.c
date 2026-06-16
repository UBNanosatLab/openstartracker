#define OST_EXPORT
#include "ost.h"

#include <png.h>

OST_DEF int ost_png_dimensions(const char *filename, int *width, int *height)
{
    FILE *fp;
    png_byte sig[8];
    png_structp png;
    png_infop info;
    png_uint_32 w;
    png_uint_32 h;

    if (!filename || !width || !height)
        return -1;
    fp = fopen(filename, "rb");
    if (!fp)
        return -2;
    if (fread(sig, 1, sizeof(sig), fp) != sizeof(sig) ||
        png_sig_cmp(sig, 0, sizeof(sig))) {
        fclose(fp);
        return -3;
    }
    png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        return -4;
    }
    info = png_create_info_struct(png);
    if (!info) {
        png_destroy_read_struct(&png, NULL, NULL);
        fclose(fp);
        return -4;
    }
    if (setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return -5;
    }
    png_init_io(png, fp);
    png_set_sig_bytes(png, (int)sizeof(sig));
    png_read_info(png, info);
    w = png_get_image_width(png, info);
    h = png_get_image_height(png, info);
    png_destroy_read_struct(&png, &info, NULL);
    fclose(fp);
    if (w > (png_uint_32)INT32_MAX || h > (png_uint_32)INT32_MAX)
        return -6;
    *width = (int)w;
    *height = (int)h;
    return 0;
}

OST_DEF int ost_png_read_rgba(const char *filename, unsigned char *dst,
                              int width, int height, int stride)
{
    FILE *fp;
    png_byte sig[8];
    png_structp png;
    png_infop info;
    png_uint_32 w;
    png_uint_32 h;
    int bit_depth;
    int color_type;
    int interlace;

    if (!filename || !dst || width <= 0 || height <= 0 || stride < 4 * width)
        return -1;
    fp = fopen(filename, "rb");
    if (!fp)
        return -2;
    if (fread(sig, 1, sizeof(sig), fp) != sizeof(sig) ||
        png_sig_cmp(sig, 0, sizeof(sig))) {
        fclose(fp);
        return -3;
    }
    png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fclose(fp);
        return -4;
    }
    info = png_create_info_struct(png);
    if (!info) {
        png_destroy_read_struct(&png, NULL, NULL);
        fclose(fp);
        return -4;
    }
    if (setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return -5;
    }

    png_init_io(png, fp);
    png_set_sig_bytes(png, (int)sizeof(sig));
    png_read_info(png, info);

    w = png_get_image_width(png, info);
    h = png_get_image_height(png, info);
    if (w != (png_uint_32)width || h != (png_uint_32)height) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return -6;
    }

    bit_depth = png_get_bit_depth(png, info);
    color_type = png_get_color_type(png, info);
    interlace = png_get_interlace_type(png, info);
    if (interlace != PNG_INTERLACE_NONE) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return -7;
    }

    if (bit_depth == 16)
        png_set_strip_16(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);
    if (png_get_valid(png, info, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png);
    if (color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);
    if (!(color_type & PNG_COLOR_MASK_ALPHA))
        png_set_filler(png, 0xff, PNG_FILLER_AFTER);

    png_read_update_info(png, info);
    if (png_get_rowbytes(png, info) != (png_size_t)(4 * width)) {
        png_destroy_read_struct(&png, &info, NULL);
        fclose(fp);
        return -8;
    }

    for (int y = 0; y < height; y++)
        png_read_row(png, (png_bytep)dst + (size_t)y * (size_t)stride, NULL);
    png_read_end(png, info);
    png_destroy_read_struct(&png, &info, NULL);
    fclose(fp);
    return 0;
}
