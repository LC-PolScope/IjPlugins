package edu.mbl.cdp.psffrequencyfilter;


import net.imglib2.Cursor;
import net.imglib2.algorithm.Benchmark;
import net.imglib2.algorithm.OutputAlgorithm;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.util.Util;


/*
 * #%L
 * ImgLib2: a general-purpose, multidimensional image processing library.
 * %%
 * Copyright (C) 2009 - 2012 Stephan Preibisch, Stephan Saalfeld, Tobias
 * Pietzsch, Albert Cardona, Barry DeZonia, Curtis Rueden, Lee Kamentsky, Larry
 * Lindsey, Johannes Schindelin, Christian Dietz, Grant Harris, Jean-Yves
 * Tinevez, Steffen Jaensch, Mark Longair, Nick Perry, and Jan Funke.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the 
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public 
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
/**
 * TODO
 *
 */
public class SmoothBandpass<T extends NumericType<T>> implements OutputAlgorithm<Img<ComplexFloatType>>, Benchmark {

	String errorMessage = "";
	boolean inPlace, bandPass;
	Img<ComplexFloatType> img, output;
	private int beginRadius, endRadius;
	private double inpix;  //e.g. = 7.4 / (100 * 1.5); // 100x objective, 1.5x optovar, 7.4um pixel of Qimaging Camera.
	private double emlambda; //e.g. = 0.598; // FM dye when bound to phospholipid membrane emits at 598nm.
	private double NA; //e.g. = 1.4; // Numerical aperture of the imaging lens.
	private double fcut; // The spatial frequency cutoff of the imaging system.
	private int sharpness;
	private int IntFac;
	//	
	private long processingTime;
	private long[] origin;

	public SmoothBandpass( //final Img<ComplexFloatType> img, final int beginRadius, final int endRadius )
			final Img<ComplexFloatType> img, double inpix, double emlambda, double NA, int sharpness, int IntFac) {

		this.img = img;
		this.inpix = inpix;
		this.emlambda = emlambda;
		this.NA = NA;
		this.sharpness = sharpness;
		this.IntFac = IntFac;
		//this.beginRadius = beginRadius;
		//this.endRadius = endRadius;
		this.inPlace = false;
		this.bandPass = true;

		fcut = 2 * NA / emlambda; // The spatial frequency cutoff of the imaging system.
		// Generate frequency grid.
		// double mcut = 1 / (2 * inpix); // Cut-off of frequency grid.
		// mx=linspace(-mcut,mcut,xlen);
		// my=linspace(-mcut,mcut,ylen);
		// [mxx, myy]=meshgrid(mx,my);
		// mrr=sqrt(mxx.^2+myy.^2);

		// double mrr=Math.sqrt(Math.pow(x,2)+Math.pow(y,2));

		// Use super-gaussian as our frequency filter.
		//FiltFreq=exp(-(mrr/fcut).^sharpness);

		//
		this.origin = new long[img.numDimensions()];
		this.origin[ 0] = img.dimension(0) - 1;
		for (int d = 1; d < this.origin.length; ++d) {
			origin[ d] = img.dimension(d) / 2;
		}
	}

	public void setImage(final Img<ComplexFloatType> img) {
		this.img = img;
	}

	public void setInPlace(final boolean inPlace) {
		this.inPlace = inPlace;
	}

	public void setBandPass(final boolean bandPass) {
		this.bandPass = bandPass;
	}

	public void setOrigin(final long[] position) {
		this.origin = position.clone();
	}

	public void setBandPassRadius(final int beginRadius, final int endRadius) {
		this.beginRadius = beginRadius;
		this.endRadius = endRadius;
	}

	public Img<ComplexFloatType> getImage() {
		return img;
	}

	public boolean getInPlace() {
		return inPlace;
	}

	public int getBeginBandPassRadius() {
		return beginRadius;
	}

	public int getEndBandPassRadius() {
		return endRadius;
	}

	public long[] getOrigin() {
		return origin;
	}

	@Override
	public boolean process() {
		final long startTime = System.currentTimeMillis();
		final Img<ComplexFloatType> img;

		if (inPlace) {
			img = this.img;
		} else {
			this.output = this.img.copy();
			img = this.output;
		}



		final long[] pos = new long[img.numDimensions()];

		final boolean actAsBandPass = bandPass;
		final Cursor<ComplexFloatType> cursor = img.cursor();
		while (cursor.hasNext()) {
			cursor.fwd();
			cursor.localize(pos);

			final float dist = Util.computeDistance(origin, pos);
			//FiltFreq=	
			double real = cursor.get().getRealDouble();
			double imaginary = cursor.get().getImaginaryDouble();
			cursor.get().setReal(          real*Math.exp(-Math.pow((dist / fcut), sharpness)));
			cursor.get().setImaginary(imaginary*Math.exp(-Math.pow((dist / fcut), sharpness)));
			
			// ++ beyond a certain radius, just set to zero ?
//			if (actAsBandPass) {
//				if (dist < beginRadius || dist > endRadius) {
//					cursor.get().setZero();
//				}
//			} else {
//				if (dist >= beginRadius && dist <= endRadius) {
//					cursor.get().setZero();
//				}
//			}
		}

		processingTime = System.currentTimeMillis() - startTime;

		// finished applying bandpass
		return true;
	}

	@Override
	public Img<ComplexFloatType> getResult() {
		if (inPlace) {
			return img;
		}
		return output;
	}

	@Override
	public long getProcessingTime() {
		return processingTime;
	}

	@Override
	public boolean checkInput() {
		return true;
	}

	@Override
	public String getErrorMessage() {
		return errorMessage;
	}

}
