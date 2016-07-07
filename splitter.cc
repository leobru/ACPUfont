#include <cstdio>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <getopt.h>

using std::vector;

int hsize = 21;
int vsize = 33;

class Picture : vector<double> {
    double & p(int r, int c) { return (*this)[r*width+c]; }
    struct column_iterator {
        int stride;
        double * ptr;
        column_iterator(int s = 0, double * p = NULL) : stride(s), ptr(p) { }
        double & operator*() { return *ptr; }
        column_iterator operator++() { ptr += stride; return *this; }
        column_iterator operator--() { ptr -= stride; return *this; }
        column_iterator operator+(int i) { return column_iterator(stride, ptr+i*stride); }
        column_iterator operator-(int i) { return column_iterator(stride, ptr-i*stride); }
        bool operator==(const column_iterator & other) const {
            return ptr == other.ptr;
        }
        bool operator!=(const column_iterator & other) const {
            return ptr != other.ptr;
        }
    };
    column_iterator column_begin(int n) { return column_iterator(width, &p(0, n)); }
    column_iterator column_end(int n) { return ++column_iterator(width, &p(height-1, n)); }
    
    FILE * runs[128];
    int currun;
    
public:
    const int width, height, full;
    const char * const postproc;

    Picture(const char * fname, const char * p) :
        currun(0),
        width(-1), height(-1), full(-1), postproc(p) {
        for (FILE ** it = runs; it != runs+128; ++it) *it = NULL;
	FILE * f = fname ? fopen(fname, "r") : stdin;
	if (!f) {
            fprintf(stderr, "Cannot open %s\n", fname);
            exit(1);
	}
	char magic[3];
	if (1 != fscanf(f, "%2s", magic) || strcmp("P5", magic)) {
            fprintf(stderr, "Want to see binary PGM (P5)\n");
            exit(1);
	}
	if (3 != fscanf(f, "%d %d %d",
                        const_cast<int*>(&width),
                        const_cast<int*>(&height),
                        const_cast<int*>(&full))
            || width <= 0 || height <= 0) {
            fprintf(stderr, "Could not parse PGM header\n");
            exit(1);
	}
	getc(f); // consume linefeed
	reserve(height*width);
	for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j) {
		int pixel = getc(f);
		if (pixel == EOF) {
                    fprintf(stderr, "File ended prematurely\n");
                    exit(1);
		}
		push_back(pixel);
            }
	fclose (f);
    }

    void dump(FILE * f) {
	fprintf(f, "P5\n%d %d\n%d\n", width, height, full);
	for (unsigned i = 0; i < size(); ++i) {
            putc(int((*this)[i]), f);
	}
        fclose(f);
    }

    void dump(int cline, int ccol, int lineno, int colno, int where) {
        if (cline-vsize/2 < 0 || ccol-hsize/2 < 0 ||
            cline + vsize/2 >= height || ccol + hsize/2 >= width) {
            fprintf(stderr, "Cell [%d,%d] with center (%d,%d) is out of bounds\n", lineno, colno, cline, ccol);
            return;
        }
        char fname[8192];
        sprintf(fname, "gost%03o", where);
        if (-1 == mkdir(fname, 0700) && errno != EEXIST) {
            perror(fname);
            exit(1);
        }
        if (postproc) 
            sprintf(fname, "%s > gost%03o/l%02dc%03d.pgm", postproc, where, lineno, colno);
        else 
            sprintf(fname, "gost%03o/l%02dc%03d.pgm", where, lineno, colno);
        FILE * f;
        if (postproc) {
            if (runs[currun]) {
                int ret = pclose(runs[currun]); 
                if (ret) {
                    fprintf(stderr, "Postproc failed with %d\n", ret);
                }
            }
            f = runs[currun] = popen(fname, "w") ;
            currun = (currun+1)%128;
        } else {
            f = fopen(fname, "w");
        }
        if (f == NULL) { perror(fname); exit(1); }
        fprintf(f, "P5\n%d %d\n%d\n", hsize, vsize, full);
        for (int i = cline - vsize/2; i <= cline + vsize/2; ++i) {
            for (int j = ccol - hsize/2; j <= ccol + hsize/2; ++j)
                putc(int(p(i,j)), f);
        }
        if (!postproc)
            fclose(f);
        else
            fflush(f);
    }

    void dump_marked(vector<int> * hor, vector<int> * ver) {
        FILE * f = fopen("marked.ppm", "w");
        fprintf(f, "P6\n%d %d\n%d\n", width, height, full);
        for (int i = 0; i < height; ++i) for (int j = 0; j < width; ++j) {
                int r, b, g;
                r = g = b = p(i,j);
                if (hor && std::find(hor->begin(), hor->end(), i) != hor->end()) { g = b = 0; r = full; }
                if (ver && std::find(ver->begin(), ver->end(), j) != ver->end()) { b = 0; r = g ? 0 : full; g = full; }
                putc(r, f); putc(g, f); putc(b, f);
            }
        fclose(f);
    }


    double border(int cline, int ccol, int width) {
        int sum = 0;
        for (int i = 0; i < vsize; ++i)
            for (int j = 0; j < hsize; ++j) {
                if (i < width || i + width >= vsize || j < width || j + width >= hsize)
                    sum += p(cline-vsize/2+i,ccol-hsize/2+j);
            }
        return sum;
    }

    double rectavg(int cline, int ccol, int hspread, int vspread) {
        double sum = 0;
        for (int i = -hspread/2; i <= hspread/2; ++i) {
            for (int j = -vspread/2; j <= vspread/2; ++j) {
                if (cline+i >= 0 && cline+i < height && ccol+j >= 0 && ccol+j < width) {
                    sum += p(cline+i,ccol+j);
                }
            }
        }
        return sum/(hspread|1)/(vspread|1);
    }


    void improve(int & cline, int & ccol) {
        const int wiggle = 3;
        // Try to wiggle the center to maximize whiteness of the wiggle-wide border
        int best = border(cline, ccol, wiggle);
        int besti = 0, bestj = 0;
        for (int i = -wiggle; i <= wiggle; ++i)
            for (int j = -wiggle; j <= wiggle; ++j) {
                if (i || j) {
                    int cur = border(cline+i, ccol+j, wiggle);
                    if (cur > best) {
                        best = cur; besti=i; bestj=j;
                    }	
                }
            }
        cline += besti;
        ccol += bestj;
    }

    vector<double> average_lines(int start, int finish) {
        vector<double> avgs(height);
        for (int i = 0; i < height; ++i) {
            avgs[i] = std::accumulate(&p(i, start), &p(i, finish+1), 0.0) / (double) (finish-start+1);
        }
        return avgs;
    }

    vector<double> average_columns(int start, int finish) {
        vector<double> avgs(width);
        for (int i = 0; i < width; ++i) {
            avgs[i] = std::accumulate(column_begin(i)+start, column_begin(i)+finish+1, 0.0) / (double) (finish-start+1);
        }
        return avgs;
    }

    ~Picture() {
        for (; runs[currun]; currun = (currun+1)%128) { pclose(runs[currun]); runs[currun] = NULL; }
    }
    
};

double smooth(const double * begin, const double * end) {
    return std::accumulate(begin, end, 0.0)/(end-begin);
}

vector<double> smooth(const vector<double> & vec, int n) {
    vector<double> smoothed(vec.size());
    for (int i = 0; i < vec.size(); ++i) {
        smoothed[i] = smooth(&vec[std::max(0, i-n)], &vec[std::min(i+n, (int)vec.size()-1)]);
    }
    return smoothed;
}
vector<int> centers(const vector<double> & vec, int cellsize) {
    vector<int> centers;
    for (int i = cellsize/2; i < vec.size()-cellsize/2; ++i) {
        if (vec[i] < vec[i-1] && vec[i] < vec[i+1]) centers.push_back(i);
    }
    return centers;
}

int main(int argc, char**argv) {
    int full, full_thr, cir_thr;
    char p2[3];
    int top_heavy_corr = 0;
    int opt;
    const char * postproc = NULL;
    bool do_dump_marked = false;
    while ((opt = getopt(argc, argv, "v:h:t:m")) != -1) {
        switch (opt) {
        case 'v':
            vsize = atoi(optarg);
            break;
        case 'h':
            hsize = atoi(optarg);
            break;
        case 't':
            top_heavy_corr = atoi(optarg);
            break;
        case 'm':
            do_dump_marked = true;
            break;
        case '?':
            fprintf(stderr, "Usage: splitter [-t N] [-v N] [-h N] [-m] [postprocessing command]\n");
            exit(1);
        }
    }
    if (optind < argc) {
        postproc = argv[optind];
        if (!*postproc) postproc = NULL;
    }
    Picture array(NULL, postproc);
    // Find local minima (darkest spots)
    vector<int> horcenters = centers(smooth(array.average_lines(0, array.width-1), vsize/4), vsize);
    vector<int> leftcenters = centers(smooth(array.average_lines(0, array.width/2), vsize/4), vsize);
    vector<int> rightcenters = centers(smooth(array.average_lines(array.width/2+1, array.width-1), vsize/4), vsize);

    double skew = 0;
    for (int i = 0; i < horcenters.size(); ++i) {
        skew += rightcenters[i]-leftcenters[i];
    }
    skew /= horcenters.size();

    // Find average interval
    double interv = 0;
    for (uint i = 1; i < horcenters.size(); ++i) interv += horcenters[i]-horcenters[i-1];
    interv /= horcenters.size()-1;
    printf("Found %d lines, average line interval = %.2f\n", horcenters.size(), interv);

    vector<int> vercenters = centers(smooth(array.average_columns(0, array.height-1), hsize/4), hsize);
    interv = 0;
    for (uint i = 1; i < vercenters.size(); ++i) interv += vercenters[i]-vercenters[i-1];
    interv /= vercenters.size()-1;
    printf("Found %d columns, average pitch = %.2f\n", vercenters.size(), interv);
    if (do_dump_marked) {
        array.dump_marked(&horcenters, &vercenters);
        exit(0);
    }
    
    for (int line = 0; line < (int)horcenters.size(); ++line) {
        // Looking from the right, find the column of the apostrophe (that's when whiteness drops)
        uint col = vercenters.size()-1;
        double avg = array.rectavg(horcenters[line], vercenters[col], hsize/2, vsize/2);
        for (; --col; ) {
            double cur = array.rectavg(horcenters[line], vercenters[col], hsize/2, vsize/2);
            if (cur < avg*0.95) break;
            avg = cur;
        }
        // 'col' is the column number of GOST 0137
        col -= 0137;
        fprintf(stderr, "Dumping line %d starting from column %d\n", line, col);
	for (int gost = 0; gost <= 0137; ++gost) {
            int hci = horcenters[line];
            int v4 = (vercenters.size()+2)/4, v2 = (vercenters.size()+1)/2;
            int curcol = gost + col;
            int hc = curcol < v4  ? hci-skew : curcol < v2 ? hci-skew/2 : curcol < v2+v4 ? hci+skew/2 : hci+skew;
            int vc = vercenters[curcol];
            // improve(array, hc, vc);
            array.dump(hc+top_heavy_corr, vc, line, col, gost);
        }
    }
}
