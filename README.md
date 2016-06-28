# ACPUfont

Extracting a drum printer font from a low quality test printout.

## Rationale

There once was a Soviet drum printer (АЦПУ-128 (ACPU-128), standing for Алфавитно-Цифровое Печатающее Устройство, AlphaNumeric Printing Unit, 128 characters wide)  with the [GOST 10859](http://en.wikipedia.org/wiki/GOST_10859) character set, excluding a few characters with codes higher than 0137, and there is a 200 dpi scan of a [one-page sample of its diagonal test](http://mailcom.com/besm6/ACPU-128.jpg).
As you can see, the paper was a little crumpled, the ribbon hasn't been re-inked for quite a while, and the printing mechanism causes a substantial ink spread at the top of the glyphs.

I've attempted to find out if there is a font similar to the ACPU-128 font. It turned out that there is none, or at least http://www.identifont.com/ doesn't know about it. After entering the salient features, even ignoring the approximation of Y with the Cyrillic У, the only font satisfying all of them is a [comic-style font](http://www.identifont.com/identify?15+.+1KT+6ZR+1L1+11+53K+19+4Y+6XA+4C+DI+1KJ+1QY+1T+9Z).

I've decided to attempt to recreate the font by enhancing the glyphs using the fact that each of them is repeated many times.

## Segmentation

The first stage is segmenting the image into individual glyphs.  It is known that there were 72 lines per page of the standard (1 ft) fanfold paper; therefore the font size is 12 pt (~33 px in 200 dpi). The pitch is ostensibly 16 = 7.5 pt = ~21 px.

First, it makes sense to extract as much as possible from the frequency domain by using the scaling feature of the djpeg tool.
At the same time the tractor holes can be removed: 
`djpeg -dct float -scale 16/8 -gray ACPU-128.jpg | pamcut -left 400 -right 6100 > crop.pgm`

The program splitter.cc (not committed yet) looks for mean centers of glyphs vertically and horizontally, accounts for some downward drift toward the right side (either a scanning artifact or a consequence of each line being struck with 4 I/O instructions per character), and for the fact that the characters are top-heavy due to dot gain. The result is [here](http://mailcom.com/besm6/acpu/tiles.zip). 

## Averaging 

The program averager.cc takes a list of files representing samples of a single character, correlates them, and outputs the correlated result. 

## Implementation details

I didn't bother with the NetPBM library; reading and writing of PGM files (the only ones I care about) is done manually.



