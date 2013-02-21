package edu.mbl.cdp.psffrequencyfilter;


/**
 * This program is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License 2 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if
 * not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 *
 * An exception is the 1D FFT implementation of Dave Hale which we use as a library, which is
 * released under the terms of the Common Public License - v1.0, which is available at
 * http://www.eclipse.org/legal/cpl-v10.html
 *
 * @author Stephan Preibisch (stephan.preibisch@gmx.de)
 */
import ij.ImageJ;
import mpicbg.util.RealSum;
import net.imglib2.Cursor;
import net.imglib2.algorithm.fft.FourierConvolution;
import net.imglib2.algorithm.fft.FourierTransform;
import net.imglib2.algorithm.fft.InverseFourierTransform;
import net.imglib2.display.ComplexImaginaryFloatConverter;
import net.imglib2.display.ComplexPhaseFloatConverter;
import net.imglib2.display.ComplexRealFloatConverter;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.io.ImgIOException;
import net.imglib2.io.ImgOpener;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

/**
 * Perform template matching by convolution in the Fourier domain
 *
 * @author Stephan Preibisch & Stephan Saalfeld
 * @author Grant Harris, Shalin Mehta
 *
 */
/*
 *  As ImageJ Plugin...
 * 
 *  + callable by macro
 *  + operate on Stacks
 * 
 */
public class PsfFrequencyFilter {

	public void test() throws ImgIOException, IncompatibleTypeException {
		// input parameters
		double inpix = 7.4 / (100 * 1.5); // 100x objective, 1.5x optovar, 7.4um pixel of Qimaging Camera.
		double emlambda = 0.598; // FM dye when bound to phospholipid membrane emits at 598nm.
		double NA = 1.4; // Numerical aperture of the imaging lens.
		double fcut = 2 * NA / emlambda; // The spatial frequency cutoff of the imaging system.
		int sharpness = 100;  // sharpness of dropoff of Gaussian
		int IntFac = 1;  // interpolation factor
		// open with ImgOpener using an ArrayImgFactory
		Img< FloatType> image = new ImgOpener() .openImg("airypattern.tif",
				new ArrayImgFactory< FloatType>(), new FloatType());

		// display image and template
		ImageJFunctions.show(image).setTitle("input");

		// compute fourier transform of the template
		final FourierTransform< FloatType, ComplexFloatType> forwardFft =
				new FourierTransform< FloatType, ComplexFloatType>(
				image, new ComplexFloatType());
		forwardFft.process();
		final Img< ComplexFloatType> spectrum = forwardFft.getResult();

		// display fft (by default in generalized log power spectrum
		ImageJFunctions.show(spectrum).setTitle("fft power spectrum");
		// display fft phase spectrum
		ImageJFunctions.show(spectrum,
				new ComplexPhaseFloatConverter< ComplexFloatType>())
				.setTitle("fft phase spectrum");
		// display fft real values
		ImageJFunctions.show(spectrum,
				new ComplexRealFloatConverter< ComplexFloatType>())
				.setTitle("fft real values");
		// display fft imaginary values
		ImageJFunctions.show(spectrum,
				new ComplexImaginaryFloatConverter< ComplexFloatType>())
				.setTitle("fft imaginary values");
		//
		// Apply 
		SmoothBandpass<ComplexFloatType> sp = new SmoothBandpass<ComplexFloatType>(spectrum,
				inpix, emlambda, NA, sharpness, IntFac);
		sp.process();
		final Img< ComplexFloatType> filtered = sp.getResult();
		//
		// compute inverse fourier transform
		final InverseFourierTransform< FloatType, ComplexFloatType> ifft =
				new InverseFourierTransform< FloatType, ComplexFloatType>(filtered, forwardFft);
		ifft.process();
		final Img< FloatType> inverse = ifft.getResult();
		ImageJFunctions.show(inverse).setTitle("Result ");

		// normalize the inverse template
		// Example6b.norm(inverse);

		// compute fourier convolution of the inverse template and the image and display it
		// ImageJFunctions.show(FourierConvolution.convolve(image, inverse)).setTitle("Convolved");
	}


	public static void main(String[] args) throws ImgIOException, IncompatibleTypeException {
		// open an ImageJ window
		new ImageJ();

		// run the example
		new PsfFrequencyFilter().test();
	}

	/**
	 * Norms all image values so that their sum is 1
	 *
	 * @param iterable - the image data
	 */
	public static void norm(final Iterable< FloatType> iterable) {
		final double sum = sumImage(iterable);

		for (final FloatType type : iterable) {
			type.setReal(type.get() / sum);
		}
	}

	/**
	 * Computes the sum of all pixels in an iterable using RealSum
	 *
	 * @param iterable - the image data
	 * @return - the sum of values
	 */
	public static < T extends RealType< T>> double sumImage(final Iterable< T> iterable) {
		final RealSum sum = new RealSum();

		for (final T type : iterable) {
			sum.add(type.getRealDouble());
		}

		return sum.getSum();
	}

	/*
	 * PsfDenoise
	 * function [out,outscale]=imFilterInterpolateMicData(in,inscale,fcut,sharpness,IntFac) 

	 % imFilterInterpolateMicData Filter and interpolate microscopy image while avoiding artifacts.
	 %
	 % Images acquired with lenses (microscope, camera) can possess fine
	 % features only upto the spatial frequency cut-off of the optics. Any fine
	 % features beyond that are due to noise. A simple and effective noise
	 % removal strategy is to remove the above-cutoff spatial frequencies.
	 % 
	 % imFilterInterpolateMicData implements filtering of the imaging data with
	 % super-gaussian filter in such a way that filtering artifacts are
	 % minimized.
	 % 
	 % The image can be optionally interpolated in frequency domain while
	 % maintaining physical variation in intensity.
	 % 
	 % USAGE: [out, outpix] =
	 % imFilterInterpolateMicData(in,inpix,fcut,sharpness,IntFac)
	 % 
	 % OUTPUTS: 
	 % out - filtered and interpolated output image. 
	 % outpix - pixel-size in the output image.
	 %
	 % INPUTS: 
	 % in - raw image 
	 % inpix - pixel-size in the input image (numerical
	 % value in your chosen units of distance). 
	 % fcut - spatial-frequency cutoff
	 % (numerical value in the inverse units of distance). 
	 % sharpness - power of the super-gaussian filter that controls transition
	 % from 1 to 0 around the cutoff. Higher value leads to faster roll-off.
	 % sharpness=2 corrsponds to the standard Gaussian filter.
	 % IntFac - Interpolation factor. 
	 % 
	 % Author and Copyright: Shalin Mehta (www.mshalin.com)
	 % License: BSD 
	 % Version history: April 2012, initial implementation with gaussian filter.
	 %                  August 2012, use super-gaussian filter.
	 %                  Feb 2013, added functionality to interpolate in
	 %                  frequency domain.
	 */
}
