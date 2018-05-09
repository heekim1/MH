package com.lifetech.converge.microhap.engine;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
//import java.util.StringJoiner;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import com.lifetech.converge.microhap.domain.HotSpot;
import com.lifetech.converge.microhap.domain.Marker;
import com.lifetech.converge.microhap.domain.OrderBy;
import com.lifetech.converge.microhap.domain.SortBy;
import com.lifetech.converge.microhap.domain.Allele;
import com.lifetech.converge.microhap.domain.Criteria;



public class Microhaplotyper {
	
	private MicrohaplotypeHotSpots hotSpots;          //hash of markerId as key and instance of HotSpot as value 
	private MicrohaplotypeMarkers markers;            //array of instance of Marker
	private SamReader samReader;                      //iterator of SAMRecord
	
	private Map<String,Map<String,Allele>> microhaplotypes;
	private int readCount;
	private int shortReadCount;
	private Criteria<Allele> criteria;
	private Float minCovThreshold;
	
	public Microhaplotyper(File hotSpotBedFile, File bamFile) throws FileNotFoundException {
		
		hotSpots = new  MicrohaplotypeHotSpots(hotSpotBedFile);
		samReader = SamReaderFactory.makeDefault().open(bamFile);
		markers = null;
		
		microhaplotypes = new HashMap<String,Map<String,Allele>>();
		readCount = 0;
		shortReadCount = 0;
		
	}
	
	
	public Microhaplotyper(File hotSpotBedFile, File bamFile, File markerFile) throws FileNotFoundException{
		
		this(hotSpotBedFile, bamFile);
		markers = new MicrohaplotypeMarkers(markerFile);
		
	}
	
	public Microhaplotyper run(){
		countAlleleCoverage();
		print();
		return this;
	}
	
	
	private void countAlleleCoverage(){
		
		for(String markerId: hotSpots.keySet()){                              //for each marker:
			
			SortedSet<HotSpot> locus = hotSpots.get(markerId);
			Map<String,Allele> alleles = new HashMap<String,Allele>();			
			SAMRecordIterator reads = getPileup(locus);                      //pileup
			
			while (reads.hasNext()) {                                        //for each read:
				SAMRecord r = reads.next();
	            String mappingSeq = getMappingSequence(r);                   //handle CIGAR string
	            try{
	            	String mhSeq = getMicrohapSeq(locus,mappingSeq,r);       //get allele (microhap)
	            	readCount++;
	            	if(alleles.containsKey(mhSeq)){                          //update allele (microhap) frequency
	            		Allele allele = alleles.get(mhSeq);
	            		allele.setalleleCoverage(allele.getAlleleCoverage()+1);
	            		allele.setTotalMQ(allele.getTotalMQ()+r.getMappingQuality());
	            		if(r.getReadNegativeStrandFlag()){
	            			allele.setMinusStrand(allele.getMinusStrand()+1);
	            		}else{
	            			allele.setPlusStrand(allele.getPlusStrand()+1);
	            		}
	            		
	            	}else{
	            		Allele allele = new Allele(mhSeq);
	            		allele.setalleleCoverage(1);
	            		allele.setTotalMQ(r.getMappingQuality());
	            		alleles.put(mhSeq, allele);
	            		if(r.getReadNegativeStrandFlag()){
	            			allele.setMinusStrand(1);
	            		}else{
	            			allele.setPlusStrand(1);
	            		}
	            	}
	            }catch(StringIndexOutOfBoundsException e){
	            	shortReadCount++;
	            }
				
	        }
			reads.close();
			updateAllelleStats(alleles);
			microhaplotypes.put(markerId, alleles);                         
		}
		
	}
	
	private void print(){
		// print header
		System.out.println("##readCount : " + readCount + ", ShortReadCount : " + shortReadCount);
		System.out.println("#MarkerId\tAvgAe\tProbOfDetMix\tallele\tcoverage\tminus_cov\tplus_cov\tavgMQ\tSNPs");
		// print body
		if(markers != null){ 
			markers.setSortBy(SortBy.Ae, OrderBy.DESC);
			for(Marker marker: markers){
				// the loci (markers) in the microhap info should contain the loci in the hostpot. 
				if(microhaplotypes.containsKey(marker.getMarkerId())){
					Map<String,Allele> alleles = microhaplotypes.get(marker.getMarkerId());
					for(String mhSeq: alleles.keySet()){
						Allele allele = alleles.get(mhSeq);
						System.out.println(marker.getMarkerId() +"\t" +  marker.getAvgAe() + "\t" + marker.getProbOfDetMix() + "\t"+
						                   allele.getMicrohap() + "\t" + allele.getAlleleCoverage() + "\t" + 
						                   allele.getMinusStrand() + "\t" + allele.getPlusStrand() + "\t" + allele.getAvgMQ() + "\t" + marker.getSNPs());
						
					}
				}
			}	
		}else{
			SortedSet<String> keys = new TreeSet<String>(microhaplotypes.keySet());
			for(String markerId: keys){
				// in this option, STDOUT will generate filtered output.
				// Map<String,Allele> alleles = (criteria != null) ? criteria.meetCriteria(microhaplotypes.get(markerId)) :
				//	                             microhaplotypes.get(markerId);
				Map<String,Allele> alleles = microhaplotypes.get(markerId);
				for(String mhSeq: alleles.keySet()){
					Allele allele = alleles.get(mhSeq);
					System.out.println(markerId +"\t"  + "\t" + "" + "\t" +
					                   allele.getMicrohap() + "\t" + allele.getAlleleCoverage() + "\t" + 
							           allele.getMinusStrand() + "\t" + allele.getPlusStrand() + "\t" + allele.getAvgMQ());
				}
			}
		}
	}
	
	public void report(File outFile) throws IOException{
		int minNumOfCnt = 0;
		int minCoverage = 0;
		
		StringBuilder sb = new StringBuilder();
		BufferedWriter br = new BufferedWriter(new FileWriter(outFile));
		
		// prepare body
		if(markers != null){  // if marker information file.
			markers.setSortBy(SortBy.Ae, OrderBy.DESC);
			for(Marker marker: markers){
				if(microhaplotypes.containsKey(marker.getMarkerId())){
					Map<String,Allele> alleles = (criteria != null) ? 
							                     criteria.meetCriteria(microhaplotypes.get(marker.getMarkerId())) :
						                         microhaplotypes.get(marker.getMarkerId());
		            int numOfCnt = getNumberOfContributors(alleles.size());
		            if(numOfCnt > minNumOfCnt){
		            	minNumOfCnt = numOfCnt;
		            }
		            
		            minCoverage = getMinCoverage(getTotalCoverage(microhaplotypes.get(marker.getMarkerId())));
		            
		            sb.append(marker.getMarkerId() + "\t" + marker.getAvgAe() + "\t" + 
			        		  marker.getProbOfDetMix() + "\t" + alleles.size() + "\t" + numOfCnt +
		                      "\t" + AllelesToString(alleles) + "\t" + AllelesToDetailString(alleles) + "\t" + marker.getSNPs() + "\t" + minCoverage + "\n");
			        

				}
	        }
		}else{              // if hotspot file due to no marker information.
			SortedSet<String> keys = new TreeSet<String>(microhaplotypes.keySet());
			for(String markerId: keys){
				Map<String,Allele> alleles = (criteria != null) ? 
						                     criteria.meetCriteria(microhaplotypes.get(markerId)) :
						                     microhaplotypes.get(markerId);
		        int numOfCnt = getNumberOfContributors(alleles.size());
		        if(numOfCnt > minNumOfCnt){
		        	minNumOfCnt = numOfCnt;
				}
		        
		        minCoverage = getMinCoverage(getTotalCoverage(microhaplotypes.get(microhaplotypes.get(markerId))));
		        
		        sb.append(markerId + "\t" + "" + "\t" + 
				         "" + "\t" + alleles.size() + "\t" + numOfCnt +
						 "\t" + AllelesToString(alleles) + "\t" + AllelesToDetailString(alleles) + "\t" + minCoverage  + "\n");
		   
			}
		}
		
		// write all
		br.write("##Filters : " + criteria.toString() + "\n");
		br.write("#The minimum number of contributors : " + minNumOfCnt + "\n");
		br.write("MarkerId\tAvgAe\tProbOfDetMix\tNumOfAlleles\tNumOfContributors\tAlleleCoverage\tAlleleCoverageDetail(-strand:+strand:avgMQ)\tSNPs\tMinimumCoverage\n");
		br.write(sb.toString());
		br.close();
	}
	
	public int getNumberOfContributors(int num){
		return (num + 2 - 1) / 2;
	}
	
	public void setFilter(Criteria<Allele> criteria){
		this.criteria = criteria;
	}
	
	public int getTotalReadCount(){
		return readCount;
	}
	
	public int getShortReadCount(){
		return shortReadCount;
	}
	
	public Float getMinCovThreshold() {
		return minCovThreshold;
	}


	public void setMinCovThreshold(Float minCovThreshold) {
		this.minCovThreshold = minCovThreshold;
	}

	private int getMinCoverage(int total){
		return (int) (total * getMinCovThreshold());
	}
	
	private int getTotalCoverage(List<Allele> alleles){
		int total = 0;
		
		for(Allele allele : alleles){
			total = total + allele.getAlleleCoverage();
		}
		return total;
	}
	
	private int getTotalCoverage(Map<String, Allele> alleles){
		int total = 0;
		
		for(String mhSeq: alleles.keySet()){
			total = total + alleles.get(mhSeq).getAlleleCoverage();
		}
		return total;
	}
	
	private String AllelesToString(Map<String,Allele> alleles){
		//StringJoiner sj = new StringJoiner(",");   //when JDK1.8
		StringBuilder sb = new StringBuilder();
		for(String allele: alleles.keySet()){
			String alleleCount = "[".concat(allele).concat("]")
					                .concat(String.valueOf(alleles.get(allele).getAlleleCoverage()));
			//sj.add(alleleCount);  //when JDK1.8
		    sb.append(alleleCount);
			sb.append(",");
		}
		//return sj.toString(); //when JDK1.8
		//return sb.substring(0, sb.length() - 1);
		return sb.length() > 0 ? sb.substring(0, sb.length() - 1) : "";
	}
	
	private String AllelesToDetailString(Map<String,Allele> alleles){
		//StringJoiner sj = new StringJoiner(",");   //when JDK1.8
		StringBuilder sb = new StringBuilder();
		for(String allele: alleles.keySet()){
			String alleleCount = "[".concat(allele).concat("]")
					                .concat(String.valueOf(alleles.get(allele).getAlleleCoverage()))
					                .concat("(").concat(String.valueOf(alleles.get(allele).getMinusStrand()))
					                .concat(":").concat(String.valueOf(alleles.get(allele).getPlusStrand()))
					                .concat(":").concat(String.valueOf(alleles.get(allele).getAvgMQ()))
					                .concat(")");
			//sj.add(alleleCount);  //when JDK1.8
		    sb.append(alleleCount);
			sb.append(",");
		}
		//return sj.toString(); //when JDK1.8
		//return sb.substring(0, sb.length() - 1);
		return sb.length() > 0 ? sb.substring(0, sb.length() - 1) : "";
	}
	
	private void updateAllelleStats(Map<String,Allele> alleles){
		updateAvgMQ(alleles);
	}
	
	private void updateAvgMQ(Map<String,Allele> alleles){
		for(String mhSeq: alleles.keySet()){
			Allele allele = alleles.get(mhSeq);
			allele.setAvgMQ(allele.getTotalMQ()/allele.getAlleleCoverage());
		}
	}
	
	private String getMicrohapSeq(SortedSet<HotSpot> locus, String seq, SAMRecord r) throws StringIndexOutOfBoundsException{
		StringBuffer sb = new StringBuffer();
		
		for(HotSpot hotspot: locus){
			int offset = hotspot.getPos() - r.getAlignmentStart();			
			sb.append(seq.charAt(offset));
		}
		
		return sb.toString();
	}
	
	
	private String getMappingSequence(SAMRecord r){
		StringBuffer readStr = new StringBuffer(r.getReadString());
		String cigarStr = r.getCigarString();
		StringBuffer newStr = new StringBuffer();
		String REGEX = "(\\d+)([SMID])"; 
		
		Pattern p = Pattern.compile(REGEX);
		Matcher m = p.matcher(cigarStr);
		while(m.find()){	
			int exclusive_stop = Integer.valueOf(m.group(1));
			if (m.group(2).equals("S") ) {
				readStr.delete(0, exclusive_stop);
			}else if (m.group(2).equals("M") ){
				newStr.append(readStr.subSequence(0, exclusive_stop));
				readStr.delete(0,exclusive_stop);
			}else if (m.group(2).equals("I") ){
				readStr.delete(0,exclusive_stop);
			}else if (m.group(2).equals("D") ){
				for(int i = 0; i < exclusive_stop; i++){
					newStr.append("_");
				}
			}
		}
		
		return newStr.toString();
	}
	
	
	public SAMRecordIterator getPileup(SortedSet<HotSpot> locus){
		
		SAMRecordIterator iterator = samReader.query(locus.first().getChr(),
				                                     locus.first().getPos(),
				                                     locus.last().getPos(),false);
		
		return iterator;
	}
}