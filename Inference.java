package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by doriswang on 5/26/17.
 */
public class Inference {
    int testCase = 20;
    int ITERATION = 1000;

    public void infer(int testruns) throws IOException, ParseException, InterruptedException{

        for(int i=0; i<testCase; i++){
            double random = Math.random()*500 - 1;
            int stNum = (int)random;

            IterativeImprovement ii = new IterativeImprovement();
            long start = System.currentTimeMillis();
            //TODO trueSeq = loadtrue data(): trueGT trueAln
            //List<Alignment> alns = ii.loadData(stNum);
            //ii.iigt(alns, ITERATION);

            long end = System.currentTimeMillis();
            long costtime = end - start;
            //TODO: print ll
            String resultFolder = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/" + ITERATION + "/" + testCase + "/" + i + "/";
            BufferedWriter llOut = new BufferedWriter(new FileWriter(resultFolder + "RunningTime.txt"));
            for (int k = 0; k < ITERATION; k++) {
                llOut.write("Running time:" + String.valueOf(costtime) + "\n");
            }
            //System.out.println("The likelihood of P(ST|S) is : " + costtime);
            System.out.println("Running time: " + costtime);
        }
    }
    public static void main(String[] args) throws IOException, ParseException, InterruptedException{
        //IterativeImprovement ii = new IterativeImprovement(trueSpecies);
        Inference in = new Inference();
        in.infer(20);
//        Network trueSpecies = new BniNetwork<NetNodeInfo>();
//        IterativeImprovement ii = new IterativeImprovement(trueSpecies);
//        long start = System.currentTimeMillis();
//
//        ii.iigt(trueSeq, ITERATION);
//
//        long end = System.currentTimeMillis();
//        long costtime = end - start;
//        //TODO: print ll
//        String resultFolder = RESULT_DIR + ITERATION + "/";
//        BufferedWriter llOut = new BufferedWriter(new FileWriter(resultFolder + "RunningTime.txt"));
//        for (int k = 0; k < ITERATION; k++) {
//            llOut.write("Running time:" + String.valueOf(costtime) + "\n");
//        }
//        //System.out.println("The likelihood of P(ST|S) is : " + costtime);
//        System.out.println("Running time: " + costtime);


    }
}
