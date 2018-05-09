package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlleleRelativeCoveragePercentCriteria implements Criteria<Allele> {
	private Float min;    //percentage
	private Float max;    //percentage

	public AlleleRelativeCoveragePercentCriteria(){
		this.min = new Float(0.0);
		this.max = new Float(1.0);
	}
	
	public List<Allele> meetCriteria(List<Allele> alleles) {
		int total = getTotalCoverage(alleles);
		int minCov = (int) (total * this.min);
		int maxCov = (int) (total * this.max);
		
		List<Allele> alleleList = new ArrayList<Allele>();
		for(Allele allele : alleles){
			if(allele.getAlleleCoverage() >= minCov && allele.getAlleleCoverage() <= maxCov){
				alleleList.add(allele);
			}
		}
		
		return alleleList;
	}
	
	public Map<String, Allele> meetCriteria(Map<String, Allele> alleles) {
		int total = getTotalCoverage(alleles);
		int minCov = (int) (total * this.min);
		int maxCov = (int) (total * this.max);
		
		Map<String,Allele> alleleMap = new HashMap<String,Allele>();
		
		for(String mhSeq: alleles.keySet()){
			Allele allele = alleles.get(mhSeq);
			if(allele.getAlleleCoverage() >= minCov && allele.getAlleleCoverage() <= maxCov){
				alleleMap.put(mhSeq, allele);
			}
		}
		return alleleMap;
	}

	public int getTotalCoverage(List<Allele> alleles){
		int total = 0;
		
		for(Allele allele : alleles){
			total = total + allele.getAlleleCoverage();
		}
		return total;
	}
	
	
	public int getTotalCoverage(Map<String, Allele> alleles){
		int total = 0;
		
		for(String mhSeq: alleles.keySet()){
			total = total + alleles.get(mhSeq).getAlleleCoverage();
		}
		return total;
	}

	
	public Float getMin() {
		return this.min;
	}

	
	public void setMin(Float min) {
		this.min = min;
	}

	
	public Float getMax() {
		return this.max;
	}

	
	public void setMax(Float max) {
		this.max = max;
	}

	
	@Override
	public String toString() {
		return "AlleleRelativeCoveragePercentCriteira [MIN=" + min + ", MAX=" + max + "]";
	}
}
