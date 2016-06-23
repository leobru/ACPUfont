/*
 * Averaging multiple unaligned grayscale pictures.
 * (c) Leo Broukhis, 2016
 */
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>

using std::vector;
using std::min;
using std::max;

class Picture : public vector<int> {
public:
    static int decontrast;
    static double msepower;
    double weight;
    int width, height;
    Picture(int w = 0, int h = 0) : vector<int>(h*w), weight(0), width(w), height(h) { }
    int & p(int i, int j) { return (*this)[i*width+j]; }
    int p(int i, int j) const { return (*this)[i*width+j]; }
    Picture(const char * fname) : weight(1) {
        FILE * f = fopen(fname, "r");
        if (!f) {
            fprintf(stderr, "Cannot open %s\n", fname);
            exit(1);
        }
        int full;
        char p2[3];
        fscanf(f, "%2s", p2);
        if (strcmp("P5", p2)) {
            fprintf(stderr, "Want to see binary PGM (P5)\n");
            exit(1);
        }
        if (3 != fscanf(f, "%d %d %d", &width, &height, &full) || full != 255) {
            fprintf(stderr, "Could not parse PGM header\n");
            exit(1);
        }
        getc(f); // eat linefeed
        reserve(height*width);
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j) {
                int pixel = getc(f);
                if (pixel == EOF) {
                    fprintf(stderr, "File ended prematurely\n");
                    exit(1);
                }
                push_back(full-pixel);
                // Now array is in degrees of blackness
            }
        fclose (f);
    }

    void dump(FILE * f) const {
	fprintf(f, "P5\n%d %d\n255\n", width, height);
	for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j)
                putc(255-(p(i,j)+weight/2)/weight, f);
	}
        fclose(f);
    }

    // Favor larger matching rectangles: the returned value is MSE divided by the number of matched pixels.
    double mse(const Picture & other, int voff, int hoff) const {
	unsigned totpix = 0;
	double sum = 0.0;
	for (int i = max(0, voff); i < min(height, height+voff); ++i) {
            for (int j = std::max(0, hoff); j < min(width, width+hoff); ++j) {
                double diff=p(i,j)/weight - other.p(i-voff,j-hoff)/other.weight;
                sum += diff*diff;
                ++totpix;
            }
        }
	return sum / pow(totpix, msepower);
    }

    Picture add(const Picture & other, int voff, int hoff) const {
	Picture ret = *this;
	for (int i = max(0, voff); i < min(height, height+voff); ++i)
	for (int j = max(0, hoff); j < min(width, width+hoff); ++j) {
            int pixel = other.p(i-voff,j-hoff);
            if (decontrast != -1) {
                // Experiment with decontrasting dot gain (useful for adding to the mean)
                double pivot = decontrast ? decontrast*other.weight : other.find_knee(true);
                if (pixel <= pivot) {
                    pixel = lround(pixel*255*other.weight/pivot);
                } else {
                    pixel = lround(255*other.weight + (pivot - pixel)/2);
                }
            }
            ret.p(i,j) += pixel; 
	}
        ret.weight += other.weight;
	return ret;
    }
    Picture divide() {
	Picture ret = *this;
	for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j)
		ret.p(i,j) = (ret.p(i,j) + (weight+1)/2)/weight;
        ret.weight = 1;
	return ret;
    }

    double best_offset(const Picture & other, int & voff, int & hoff) const {
	// Compute average distances for all offsets up to size/3
	double minmse = 10e38;
	for (int i = -height/3; i <= height/3; ++i) 
	for (int j = -width/3; j <= width/3; ++j) {
		double curmse = mse(other, i, j);
		if (curmse < minmse) {
			// fprintf(stderr, "MSE = %.3f voff=%d hoff=%d\n", curmse, i, j);
			minmse = curmse;
			voff = i;
			hoff = j;
		}
	}
	// fprintf(stderr, "Best offsets are %d %d\n", voff, hoff);
        return minmse;
    }

    Picture contrast(int pivot, int power) {
        Picture ret = *this;
        if (weight == 0.0) return ret;
        double knee = pivot/255.0;
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j) {
                double v = p(i,j)/(weight*255.0);
                if (v <= knee) {
                    v = pow(v/knee/2, power)/pow(0.5, power-1);
                } else {
                    v = 1-pow((1-v)/(1-knee)/2, power)/pow(0.5, power-1);
                }
                ret.p(i,j) = lround(v*255*weight);
            }
        return ret;
    }

// A pixel darker than threshold will be made progressively darker
// IF it is already darker than the median of its neighborhood,
// otherwise, lighter if it is already lighter.
    Picture contrast2(int threshold) {
        Picture ret = *this;
        for (int i = 0; i < height; ++i)
            for (int j = 0; j < width; ++j) {
                int pixel = p(i,j);
                if (pixel <= threshold) continue;
                vector<int> nbrs(9);
                int pos = 0;
                for (int a = -1; a <= 1; ++a) for (int b = -1; b <= 1; ++b) {
                        nbrs[pos++] =
                            p(max(0, min(height-1, i+a)), max(0, min(width-1, j+b)));
                    }
                std::sort(nbrs.begin(), nbrs.end());
                if (pixel == nbrs[8]) ret.p(i,j) = 255-(255-pixel)/2;
                if (pixel == nbrs[7]) ret.p(i,j) = 255-(255-pixel)/1.6;
                if (pixel == nbrs[6]) ret.p(i,j) = 255-(255-pixel)/1.4;
                if (pixel == nbrs[5]) ret.p(i,j) = 255-(255-pixel)/1.2;
                if (pixel == nbrs[3]) ret.p(i,j) = pixel/1.2;
                if (pixel == nbrs[2]) ret.p(i,j) = pixel/1.4;
                if (pixel == nbrs[1]) ret.p(i,j) = pixel/1.6;
                if (pixel == nbrs[0]) ret.p(i,j) = pixel/2;
            }
        return ret;
    }

    int find_knee(bool upper) const {
        vector<int> pixels = *this;
        if (upper) {
            std::sort(pixels.begin(), pixels.end(), std::greater<int>());
        } else {
            std::sort(pixels.begin(), pixels.end());
        }
        // The background/flush color ends where the normalized histogram derivative
        // is greater than 1 for the first time. 
        int totpix = pixels.size();
        int range = std::abs(pixels.back()-pixels.front());
        int grouplen = totpix/(range+1);
        // We've found the knee when the mean of the current group differs
        // from the mean of the previous group by more than weight.
        double prev = std::accumulate(pixels.begin(), pixels.begin()+grouplen, 0)/(double)grouplen;
        for (int i = 1; i < range; ++i) {
            double cur = std::accumulate(pixels.begin()+i*grouplen, pixels.begin()+(i+1)*grouplen, 0)/(double)grouplen;
            if (upper ? cur < prev-weight : cur > prev+weight) {
                 return prev;
            }
            cur = prev;
        }
        // Nothing interesting found, return the median
        return pixels[totpix/2];
    }
}; // class Picture

int Picture::decontrast = -1;
double Picture::msepower = 2;

struct offset {
    int what, with;
    int voff, hoff;
    offset(int wh, int wi, int v, int h) : what(wh), with(wi), voff(v), hoff(h) { }
    bool operator<(const offset & other ) const { return voff + hoff < other.voff + other.hoff; }
};

void averageToMean(vector<Picture> & pics, bool contrastMean) {
    Picture mean = pics.front();
    for (uint i = 1; i < pics.size(); ++i) {
        mean = mean.add(pics[i], 0, 0);
    }
    mean = mean.divide();
    mean.dump(fopen("mean.pgm", "w"));
    if (contrastMean) {
        int pivot = (mean.find_knee(false) + mean.find_knee(true))/2;
        mean = mean.contrast(pivot, 2);
        mean.dump(fopen("contrasted.pgm", "w"));
    }
    Picture & res = mean;
    // Sort by closeness to the mean
    vector<std::pair<double, int> > order(pics.size());
    for (uint i = 0; i < pics.size(); ++i) {
        order[i] = std::make_pair(mean.mse(pics[i], 0, 0), i);
    }
    std::sort(order.begin(), order.end());
    for (uint i = 0; i < pics.size(); ++i) { 
        int voff = 0, hoff = 0;
        res.best_offset(pics[order[i].second], voff, hoff);
        res = res.add(pics[order[i].second], voff, hoff);
    }
    res.dump(stdout);
}

void averageSpanningTree(vector<Picture> & pics) {
    Picture mean = pics.front();
    int remaining = pics.size();
    for (uint i = 1; i < pics.size(); ++i) {
        mean = mean.add(pics[i], 0, 0);
    }
    mean = mean.divide();
    vector<std::pair<double, offset> > offsets;
    // Try merging following the min spanning tree of pictures.
    for (uint i = 0; i < pics.size()-1; ++i) {
        for (uint j = i+1; j < pics.size(); ++j) {
            int voff = 0, hoff = 0;
            double best = pics[i].best_offset(pics[j], voff, hoff);
            offsets.push_back(std::make_pair(-best, offset(i, j, voff, hoff)));
        }
    }
    // Sort by negated best; closest pictures are at the end of the array.
    std::sort(offsets.begin(), offsets.end());
    while (!offsets.empty()) {
        offset best = offsets.back().second;
        if (pics[best.what].weight == 0 || pics[best.with].weight == 0) {
            offsets.pop_back();
            continue;
        }
        Picture next = pics[best.what].add(pics[best.with], best.voff, best.hoff);
        pics[best.what] = Picture();
        pics[best.with] = Picture();
        for (uint i = 0; i < pics.size(); ++i) {
            int voff, hoff;
            if (pics[i].weight) {
                double best = next.best_offset(pics[i], voff, hoff);
                offsets.push_back(std::make_pair(-best, offset(pics.size(), i, voff, hoff)));
            }
        }
        pics.push_back(next);
        --remaining;
    }
    // Re-center the result by correlating to the mean
    int voff = 0, hoff = 0;
    mean.best_offset(pics.back(), voff, hoff);
    Picture(mean.width, mean.height).add(pics.back(), voff, hoff).dump(stdout);
}

int main(int argc, char**argv) {
    // Take at least 2 files, mean mode (default) or spanning tree mode.
    bool spanning = false;
    bool contrastMean = false;
    int opt;
    static const char usage[] =
        "Usage: %s [-s] [-t] [-mF] [-dN] file1 file2 ...\n\
	-s - use spanning tree combining (default - adding to the mean)\n\
	-c - contrast the mean in the default mode\n\
	-mF - use power F of the pixel count when computing MSE, default 2.0\n\
	-dN - decontrast while adding samples, N=threshold, 0: automatic selection, nondigit: 128\n";

    while ((opt = getopt(argc, argv, "scd:m:")) != -1) {
        switch (opt) {
        case 's':
            spanning = true;
            break;
        case 'c':
            contrastMean = true;
            break;
        case 'd':
            Picture::decontrast = isdigit(*optarg) ? atoi(optarg) : 128;
            break;
        case 'm':
            Picture::msepower = atof(optarg);
            if (Picture::msepower != 0.0)
                break;
        default: /* '?' */
            fprintf(stderr, usage, argv[0]);
            exit(1);
        }
    }
    if (argc - optind < 2) {
        fprintf(stderr, usage, argv[0]);
        exit(1);
    }

    vector<Picture> pics;
    for (int arg = optind; arg < argc; ++arg) {
        pics.push_back(Picture(argv[arg]));
        if (pics.back().height != pics.front().height || pics.back().width != pics.front().width) {
            fprintf(stderr, "Files must have same dimensions, got %dx%d and %dx%d\n",
                    pics.front().width, pics.front().height, pics.back().width, pics.back().height);
            exit(1);
        }
    }
    
    if (spanning) {
        averageSpanningTree(pics);
    } else {
    averageToMean(pics, contrastMean);
    }
}
