package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class AlleleCoverageCriteria implements Criteria<Allele> {
	
	private int min;
	private int max;

	public AlleleCoverageCriteria(){
		this.min = 0;
		this.max = 100000;
	}
	
	public List<Allele> meetCriteria(List<Allele> alleles) {
		List<Allele> alleleList = new ArrayList<Allele>();
		for(Allele allele : alleles){
			if(allele.getAlleleCoverage() >= min && allele.getAlleleCoverage() <= max){
				alleleList.add(allele);
			}
		}
		
		return alleleList;
	}
	
	public Map<String, Allele> meetCriteria(Map<String, Allele> alleles) {
		Map<String,Allele> alleleMap = new HashMap<String,Allele>();
		
		for(String mhSeq: alleles.keySet()){
			Allele allele = alleles.get(mhSeq);
			if(allele.getAlleleCoverage() >= min && allele.getAlleleCoverage() <= max){
				alleleMap.put(mhSeq, allele);
			}
		}
		return alleleMap;
	}

	
	public int getMin() {
		return this.min;
	}

	public void setMin(int min) {
		this.min = min;
	}

	public int getMax() {
		return this.max;
	}

	public void setMax(int max) {
		this.max = max;
	}

	@Override
	public String toString() {
		return "AlleleCoverageCriteria [MIN=" + min + ", MAX=" + max + "]";
	}
}
