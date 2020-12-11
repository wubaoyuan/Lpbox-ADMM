from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef("""
    ssize_t read (int __fd, void *__buf, size_t __nbytes);
    ssize_t write (int __fd, const void *__buf, size_t __n);
    static void sendTriplet(int i, int j, double v, int fd);
""")

ffibuilder.set_source("_CInterface", r"""
    #include <unistd.h>
    struct Trip {
        int row;
        int col;
        double val;
    };
    static void sendTriplet(int i, int j, double v, int fd) {
        struct Trip triplet;
        triplet.row = i;
        triplet.col = j;
        triplet.val = v;
        write(fd, (void *)&triplet, sizeof(struct Trip));
    }
""")

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)