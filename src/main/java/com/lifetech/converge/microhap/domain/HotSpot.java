package com.lifetech.converge.microhap.domain;


public class HotSpot implements Comparable<HotSpot>{

	private String rsId;
	private final String chr;
	private final int pos;
	
	
	public HotSpot(String rsId, String chr, int pos){
		this.rsId = rsId;
		this.chr = chr;
		this.pos = pos;
	}

	
	
	public String getRsId() {
		return rsId;
	}



	public void setRsId(String rsId) {
		this.rsId = rsId;
	}



	public String getChr() {
		return chr;
	}


	public int getPos() {
		return pos;
	}

	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + pos;
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		HotSpot other = (HotSpot) obj;
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (pos != other.pos)
			return false;
		return true;
	}


	public int compareTo(HotSpot o) {
		int pos1 = pos;
		int pos2 = o.pos;
		
		if(pos1 > pos2){
			return 1;
		}else if(pos1 < pos2){
			return -1;
		}else{
			return 0;
		}
			
	}
}
