# Configuration variables governing the build of gretl for win32

# directory containing the cross tools
CROSSDIR = /opt/cross-tools/mingw32

# prefix for mingw tools (e.g. mingw32-gcc)
MGW_PREFIX = mingw32

# mingw include dir
MGW_INC = $(CROSSDIR)/include

# glib includes: adjust to match your system
GLIB_INC = $(CROSSDIR)/include/glib-2.0

# msgfmt command for producing win32 messages file
WIN32_MSGFMT = wine c:/bin/msgfmt.exe

# pkgconfig path
PKG_CONFIG_PATH = $(CROSSDIR)/lib/pkgconfig
