package com.lifetech.converge.microhap.domain;

public class Allele{
	private String microhap;
	private int alleleCoverage;
	private int totalMQ;
	private int plusStrand;
	private int minusStrand;
	private int avgMQ;
	private float StrandBias;
	
	public Allele(String microhap){
		this.microhap = microhap;
		this.alleleCoverage = 0;
		this.totalMQ = 0;
		this.plusStrand = 0;
		this.minusStrand = 0;
		this.avgMQ = 0;
		this.StrandBias = 0;
		
	}

	public String getMicrohap() {
		return microhap;
	}

	public void setMicrohap(String microhap) {
		this.microhap = microhap;
	}

	public int getAlleleCoverage() {
		return alleleCoverage;
	}

	public void setalleleCoverage(int alleleCoverage) {
		this.alleleCoverage = alleleCoverage;
	}

	public int getTotalMQ() {
		return totalMQ;
	}

	public void setTotalMQ(int totalMQ) {
		this.totalMQ = totalMQ;
	}

	public int getPlusStrand() {
		return plusStrand;
	}

	public void setPlusStrand(int plusStrand) {
		this.plusStrand = plusStrand;
	}

	public int getMinusStrand() {
		return minusStrand;
	}

	public void setMinusStrand(int minusStrand) {
		this.minusStrand = minusStrand;
	}

	public int getAvgMQ() {
		return avgMQ;
	}

	public void setAvgMQ(int avgMQ) {
		this.avgMQ = avgMQ;
	}

	public float getStrandBias() {
		return StrandBias;
	}

	public void setStrandBias(float StrandBias) {
		this.StrandBias = StrandBias;
	}

	
}