package com.lifetech.converge.microhap.domain;

public enum OrderBy {
	ASC("ascending"), DESC("descending");
	private String name;
	
	OrderBy(String name){
		this.name= name;
	}
	
	public String getName() {
        return name;
    }
}
