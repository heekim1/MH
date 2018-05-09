package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlleleSeqContainDeletionCriteria implements Criteria<Allele> {
	
	private String del;

	public AlleleSeqContainDeletionCriteria(){
		this.del = "_";
	}
	
	public List<Allele> meetCriteria(List<Allele> alleles) {
		List<Allele> alleleList = new ArrayList<Allele>();
		for(Allele allele : alleles){
			if(!allele.getMicrohap().contains(del)){
				alleleList.add(allele);
			}
		}
		
		return alleleList;
	}
	
	public Map<String, Allele> meetCriteria(Map<String, Allele> alleles) {
		Map<String,Allele> alleleMap = new HashMap<String,Allele>();
		
		for(String mhSeq: alleles.keySet()){
			Allele allele = alleles.get(mhSeq);
			if(!allele.getMicrohap().contains(del)){
				alleleMap.put(mhSeq, allele);
			}
		}
		return alleleMap;
	}

	@Override
	public String toString() {
		return "AlleleSeqContainDeletionCriteria [DEL=" + del + "]";
	}
}
