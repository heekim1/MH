package com.lifetech.converge.microhap.domain;


import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * microhaplotype marker information
 */

public class Marker {
	private final String markerName;
	private final String markerId;
	private final float probOfDetMix;
	private final int detMixRank;
	private final float avgAe;
	private final int AeRank;
	private final String chr;
	private final int pMostPos;
	
	private final String SNPs;
	
	public static class Builder {
		private String markerName;
		private String markerId;
		private float probOfDetMix;
		private int detMixRank;
		private float avgAe;
		private int AeRank;
		private String chr;
		private int pMostPos;
		private String snps;
		
		public Builder markerName(String val){
			markerName = val;
			return this;
		}
		
		public Builder markerId(String val){
			markerId = val;
			return this;
		}
		
		
		public Builder probOfDetMix(float val){
			probOfDetMix = val;
			return this;
		}
		
		public Builder detMixRank(int val){
			detMixRank = val;
			return this;
		}
		
		public Builder avgAe(float val){
			avgAe = val;
			return this;
		}
		
		public Builder AeRank(int val){
			AeRank = val;
			return this;
		}
		
		public Builder chr(String val){
			chr = val;
			return this;
		}
		
		public Builder pMostPos(int val){
			pMostPos = val;
			return this;
		}
		
		public Builder snps(String val){
			snps = val;
			return this;
		}
		
		public Marker build() {
			return new Marker(this);
		}
	}
	
	
	private Marker(Builder builder){
		markerName = builder.markerName;
		markerId = builder.markerId;
		probOfDetMix = builder.probOfDetMix;
		detMixRank = builder.detMixRank;
		avgAe = builder.avgAe;
		AeRank = builder.AeRank;
		chr = builder.chr;
		pMostPos = builder.pMostPos;
		
		//parse SNPs
		SNPs = builder.snps;
		
	}


	public String getMarkerName() {
		return markerName;
	}


	public String getMarkerId() {
		return markerId;
	}


	public float getProbOfDetMix() {
		return probOfDetMix;
	}


	public int getDetMixRank() {
		return detMixRank;
	}


	public float getAvgAe() {
		return avgAe;
	}


	public int getAeRank() {
		return AeRank;
	}


	public String getChr() {
		return chr;
	}


	public int getPMostPos() {
		return pMostPos;
	}


	public String getSNPs() {
		return SNPs;
	}
}