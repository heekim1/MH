package com.lifetech.converge.microhap.engine;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.lifetech.converge.microhap.domain.Allele;
import com.lifetech.converge.microhap.domain.AlleleAvgMQCriteria;
import com.lifetech.converge.microhap.domain.AlleleCoverageCriteria;
import com.lifetech.converge.microhap.domain.AlleleRelativeCoveragePercentCriteria;
import com.lifetech.converge.microhap.domain.AlleleSeqContainDeletionCriteria;
import com.lifetech.converge.microhap.domain.AndCriteria;
import com.lifetech.converge.microhap.domain.Criteria;

import org.apache.commons.cli.DefaultParser;


/**
 * Hello world!
 *
 */

public class Main 
{
    public static void main( String[] args ) throws ParseException, IOException
    {	
    	Options options = getOptions();
    	CommandLineParser parser = new DefaultParser();
    	CommandLine cmd = parser.parse(options, args);
    	if(checkCommand(cmd)){
    		
    	}

    	File inputBed = new File(cmd.getOptionValue("bed"));
    	File inputBam = new File(cmd.getOptionValue("bam"));
    	File inputInfo = cmd.getOptionValue("info") != null  ? new File(cmd.getOptionValue("info")) : null;
    	File outputFile = new File(cmd.getOptionValue("o"));
        
        // instantiate Microhaplotyper
        Microhaplotyper mher = (inputInfo == null) ? 
        		new  Microhaplotyper(inputBed,inputBam) : new Microhaplotyper(inputBed,inputBam, inputInfo);
        
        // define filters
        // AlleleCoverageCriteria alleleCoverageCriteria = new AlleleCoverageCriteria();
        AlleleRelativeCoveragePercentCriteria alleleRelativeCoveragePercentCriteria = new AlleleRelativeCoveragePercentCriteria();         
        AlleleAvgMQCriteria alleleAvgMQCriteria = new AlleleAvgMQCriteria();
        AlleleSeqContainDeletionCriteria alleleSeqContainDeletionCriteria= new AlleleSeqContainDeletionCriteria();
        AlleleCoverageCriteria alleleCoverageCriteria = new AlleleCoverageCriteria();        
                
        if(cmd.getOptionValue("minCov") != null)  { 
        	alleleRelativeCoveragePercentCriteria.setMin(Float.valueOf(cmd.getOptionValue("minCov")));
        	mher.setMinCovThreshold(Float.valueOf(cmd.getOptionValue("minCov")));
        }else{
        	alleleRelativeCoveragePercentCriteria.setMin(new Float(0.02));
        	mher.setMinCovThreshold(new Float(0.02));
        }
        if(cmd.getOptionValue("maxCov") != null) 
        	alleleRelativeCoveragePercentCriteria.setMax(Float.valueOf(cmd.getOptionValue("maxCov")));
       
        alleleAvgMQCriteria.setMIN(60);    // it is hard filter
        alleleCoverageCriteria.setMin(10); // it is hard filter
        
                
        Criteria<Allele> andCriteria = new AndCriteria(alleleRelativeCoveragePercentCriteria,alleleAvgMQCriteria)
                		                           .And(alleleSeqContainDeletionCriteria)
        										   .And(alleleCoverageCriteria);
        //OrCriteria orCriteria = new OrCriteria(alleleCoverageCriteria,alleleAvgMQCriteria);
        
        mher.setFilter(andCriteria);
        mher.run().report(outputFile);
        
        
        
    }
    
    public static boolean checkCommand(CommandLine cmd){
		
    	if(cmd.hasOption("b")){
    		System.out.println("b option : " + cmd.getOptionValue("b"));
    	}
    	
    	return true;
    	
    }
    
    public static Options getOptions(){
    	Option bedFile = Option.builder("bed")
    			               .hasArg()
    			               .required(true)
    			               .longOpt("bedfile")
    			               .build();
    	
    	Option bamFile = Option.builder("bam")
	               .hasArg()
	               .required(true)
	               .longOpt("bamfile")
	               .build();
    	
    	Option infoFile = Option.builder("info")
	               .hasArg()
	               .longOpt("infofile")
	               .build();
    	
    	Option outputFile = Option.builder("o")
	               .hasArg()
	               .required(true)
	               .longOpt("out")
	               .build();
    	
    	Option sId = Option.builder("s")
	               .hasArg()
	               .longOpt("sid")
	               .build();
    	
    	Option barCode = Option.builder("b")
	               .hasArg() 
	               .longOpt("barcode")
	               .build();
    	
    	Option minCovCutoff = Option.builder("minCov")
    			   .hasArg()
    			   .longOpt("minCoverageCutoff")
    			   .build();
    	
    	Option maxCovCutoff = Option.builder("maxCov")
 			   .hasArg()
 			   .longOpt("maxCoverageCutoff")
 			   .build();
    	
    	Option format = Option.builder("f")
    			   .argName("JSON|TSV|CSV")
	               .hasArg()
	               .longOpt("format")
	               .build();
    	
    	Options options = new Options();
    	options.addOption(bedFile);
    	options.addOption(bamFile);
    	options.addOption(infoFile);
    	options.addOption(sId);
    	options.addOption(barCode);
    	options.addOption(minCovCutoff);
    	options.addOption(maxCovCutoff);
    	options.addOption(outputFile);
    	options.addOption(format);
    	return options;
    			
    }
}
