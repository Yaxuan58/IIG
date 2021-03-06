package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.TemporalConstraints;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeWithBranchLengthProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Bipartitions;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import sun.nio.ch.Net;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by doriswang on 3/26/17.
 */
public class InferOperator {


    private  String InputPath;
    private  String OutputPath;
    private  Integer Index;
    private  double HalfTheta;
    private  String Iteration;
    private  String LociNum;

    protected  String _msdir;
    //    _msdir = inputPath + "tools/msFiles/msdir/ms"; ///home/yw58/IIG/tools/Astral/astral.4.11.1.jar
//    _ASTRALdir = inputPath + "tools/Astral/";
//    RESULT_DIR = outputPath ;//"/Users/doriswang/PhyloNet/Data/IIG/result/" | home/yw58/IIG/output
    protected  String _ASTRALdir;
    protected  String _RAxMLdir;
    private  String _fastTree;
    //String inputPath, String outputPath, int index, double hTheta, int iteration, int lnum
    InferOperator(String inputPath, String outputPath, int index, int rNum) {
        InputPath = inputPath;
        OutputPath = outputPath;
        _msdir = InputPath + "tools/msFiles/msdir/ms";
        //    _msdir = inputPath + "tools/msFiles/msdir/ms"; ///home/yw58/IIG/tools/Astral/astral.4.11.1.jar
//    _ASTRALdir = inputPath + "tools/Astral/";
//    RESULT_DIR = outputPath ;//"/Users/doriswang/PhyloNet/Data/IIG/result/" | home/yw58/IIG/output
        _ASTRALdir = InputPath + "tools/Astral/" + rNum + "/";
         _RAxMLdir = InputPath + "tools/RAxML/" + rNum + "/" ;
        _fastTree = InputPath + "tools/FT/"+ rNum + "/";
        Index = index;

    }

    //TODO
    public static Map<String, String> readAlignmentFromNexusFile(String inFile, int len) {
        Map<String, String> aln = new HashMap<String, String>();

        return aln;
    }

    //TODO
    public static boolean writeAlignmentNexus(Map<String, String> aln, String outFile, int firstSite, int len) {
        boolean success = false;

        return success;
    }
    public static void main(String[] args) throws  IOException, InterruptedException, ParseException {

        InferOperator ifo = new InferOperator("/Users/doriswang/PhyloNet/", "/Users/doriswang/PhyloNet/output", 0, 0);

        //Users/doriswang/Desktop/output/001/10/50/5/1000/50
        //Users/doriswang/Desktop/output/001/2/50/1/1000/GlobalBestTrees.txt
        //        ifo.MCMC_Compare();
        //ifo.gtOutfitChecker();


        //ifo.gtOutfitChecker("/Users/doriswang/Treefix/venv/treefix-1.1.10/s16/symmetricDist/200_0001_TFTrees.txt");
        Tree t1 = Trees.readTree("((((3:1.238,4:1.238):0.84,(1:1.54,2:1.54):0.54):1.718,(6:3.786,(5:2.244,(7:1.596,8:1.596):0.648):1.542):0.01):2.294,((13:1.172,14:1.172):3.072,(((11:1.498,12:1.498):0.564,(9:1.538,10:1.538):0.524):1.084,(15:2.448,16:2.448):0.7):1.096):1.848);");
        Tree t2 = Trees.readTree("((((9:1.336,10:1.336):0.866,(11:1.974,12:1.974):0.228):0.996,((15:1.162,16:1.162):1.758,(13:1.488,14:1.488):1.432):0.278):2.186,(1:4.204,(2:3.258,((3:1.222,4:1.222):1.784,(6:2.53,(5:2.046,(7:1.172,8:1.172):0.874):0.482):0.476):0.252):0.948):1.178);");
                //("((((3:1.238,4:1.238):0.84,(1:1.54,2:1.54):0.54):1.718,(6:3.786,(5:2.244,(7:1.596,8:1.596):0.648):1.542):0.01):2.294,((13:1.172,14:1.172):3.072,(((11:1.498,12:1.498):0.564,(9:1.538,10:1.538):0.524):1.084,(15:2.448,16:2.448):0.7):1.096):1.848);");
        //((((9:1.336,10:1.336):0.866,(11:1.974,12:1.974):0.228):0.996,((15:1.162,16:1.162):1.758,(13:1.488,14:1.488):1.432):0.278):2.186,(1:4.204,(2:3.258,((3:1.222,4:1.222):1.784,(6:2.53,(5:2.046,(7:1.172,8:1.172):0.874):0.482):0.476):0.252):0.948):1.178);
        System.out.print(ifo.equalTree(t1,t2));



        //ifo.gtOutfitChecker("/Users/doriswang/Desktop/output/0001/10/50/", ifo.InputPath + "input/0001/", 10, 10, 200);
        ///Users/doriswang/Treefix/venv/treefix-1.1.10/s16/a-AsmtTrees/1000_001_TFTrees.txt

    }

    // For Treefix: Check the RF_distance (trueST, TF_AST_ST)
    // Paras: UGARInfer ui = new UGARInfer("/Users/doriswang/PhyloNet/","/Users/doriswang/PhyloNet/Data/IIG/result/",1,0.005,50,5,0, 1000, 50);
    public void gtOutfitChecker() throws IOException, ParseException {
        int lociNum = 10;
        int fileNum = 20;
        Tree ast;
        String filePath = "/Users/doriswang/Treefix/venv/treefix-1.1.10/s16/";
        double[] d = new double[fileNum + 1];
        Tree trueST = Trees.readTree("(1:1.0,(2:1.0,(3:1.0,(4:1.0,(5:1.0,(6:1.0,(7:1.0,(8:1.0,(9:1.0,(10:1.0,(11:1.0,(12:1.0,(13:1.0,(14:1.0,(15:1.0,16:1.0)I14:1.0)I13:1.0)I12:1.0)I11:1.0)I10:1.0)I9:1.0)I8:1.0)I7:1.0)I6:1.0)I5:1.0)I4:1.0)I3:1.0)I2:1.0)I1:1.0)I0;");
        //(1:1.0,(2:1.0,(3:1.0,(4:1.0,(5:1.0,(6:1.0,(7:1.0,(8:1.0,(9:1.0,(10:1.0,(11:1.0,(12:1.0,(13:1.0,(14:1.0,(15:1.0,16:1.0)I14:1.0)I13:1.0)I12:1.0)I11:1.0)I10:1.0)I9:1.0)I8:1.0)I7:1.0)I6:1.0)I5:1.0)I4:1.0)I3:1.0)I2:1.0)I1:1.0)I0;

        //Tree trueST = Trees.readTree("((((16:1.0,15:1.0):1.0,(14:1.0,13:1.0):1.0):1.0,((12:1.0,11:1.0):1.0,(10:1.0,9:1.0):1.0):1.0):1.0,(((8:1.0,7:1.0):1.0,(6:1.0,5:1.0):1.0):1.0,((4:1.0,3:1.0):1.0,(2:1.0,1:1.0):1.0):1.0):1.0);");
        //Tree trueST = Trees.readTree("(((((1:1.0,2:1.0)I17:1.0,(3:1.0,4:1.0)I18:1.0)I9:1.0,((5:1.0,6:1.0)I19:1.0,(7:1.0,8:1.0)I20:1.0)I10:1.0)I5:1.0,(((9:1.0,10:1.0)I21:1.0,(11:1.0,12:1.0)I22:1.0)I11:1.0,((13:1.0,14:1.0)I23:1.0,(15:1.0,16:1.0)I24:1.0)I12:1.0)I6:1.0)I3:1.0,((((17:1.0,18:1.0)I25:1.0,(19:1.0,20:1.0)I26:1.0)I13:1.0,((21:1.0,22:1.0)I27:1.0,(23:1.0,24:1.0)I28:1.0)I14:1.0)I7:1.0,(((25:1.0,26:1.0)I29:1.0,(27:1.0,28:1.0)I30:1.0)I15:1.0,((29:1.0,30:1.0)I31:1.0,(31:1.0,32:1.0)I32:1.0)I16:1.0)I8:1.0)I4:1.0)I0;");
        double temp = 0.0;
        for (int i = 0; i < fileNum; i++) {
            List<Tree> tfTrees = new ArrayList<Tree>();
            for (int ln = 0; ln < lociNum; ln++) {
                int fileIndex = i * lociNum + ln;
                String tree = Treefix.readToString(filePath + fileIndex + "/0.nt.raxml.treefix.tree");
                Tree t = Trees.readTree(tree);
                tfTrees.add(t);
            }
            ast = Trees.readTree(initAST(tfTrees));
            d[i] = this.getDistance(ast, trueST);

            temp += d[i];
        }
        for (int i = 0; i < fileNum; i++) {
            System.out.println(d[i]);
        }
        d[fileNum] = temp / fileNum;
        System.out.println(d[fileNum]);
    }

    // For Treefix: Check the RF_distance (trueST, TF_AST_ST)
    // Paras: "/Users/doriswang/Treefix/venv/treefix-1.1.10/s16/1000_0001_TFTrees.txt "

    public void gtOutfitChecker(String tfTreesFile) throws IOException, ParseException {
        int lociNum = 10;
        int fileNum = 10;
        Tree ast;
        String filePath = tfTreesFile;
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String[] trees = reader.readLine().trim().split(";");
        double[] d = new double[fileNum + 1];
        Tree trueST = Trees.readTree("((((16:1.0,15:1.0):1.0,(14:1.0,13:1.0):1.0):1.0,((12:1.0,11:1.0):1.0,(10:1.0,9:1.0):1.0):1.0):1.0,(((8:1.0,7:1.0):1.0,(6:1.0,5:1.0):1.0):1.0,((4:1.0,3:1.0):1.0,(2:1.0,1:1.0):1.0):1.0):1.0);");
        double temp = 0.0;
        int counter = 0;
        for (int i = 0; i < fileNum; i++) {
            List<Tree> tfTrees = new ArrayList<Tree>();
            for (int ln = 0; ln < lociNum; ln++) {
                int fileIndex = i * lociNum + ln;
                String tree = trees[counter] + ";";
                Tree t = Trees.readTree(tree);
                tfTrees.add(t);
                counter++;
            }
            ast = Trees.readTree(initAST(tfTrees));
            d[i] = this.getDistance(ast, trueST);

            temp += d[i];
        }
        for (int i = 0; i < fileNum; i++) {
            System.out.println(d[i]);
        }
        d[fileNum] = temp / fileNum;
        System.out.println(d[fileNum]);
    }

    public String initAST(List<Tree> gts) throws IOException, ParseException {
//        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i /Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i " + _ASTRALdir + "test_data/tempIn.tre -o " + _ASTRALdir + "test_data/testOut.tre"; ///Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        Tree temp = gts.get(0);

        //TODO: add weight
        for (int i = 0; i < gts.size(); i++) {
            in.write(gts.get(i).toString() + "\n");
        }
        in.flush();
        in.close();
        ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
        pb.redirectErrorStream(true);
        try {
            Process proc = pb.start();
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.out.println(ex.getMessage());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        String outputFile = _ASTRALdir + "test_data/testOut.tre";
        BufferedReader out = new BufferedReader(new FileReader(outputFile));
        String stString = out.readLine().trim();
        int index0 = 0;
        int index1 = 1;
        List<Integer> indexes = new ArrayList<Integer>();
        for (int i = 1; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
            index1++;
        }
        for (int i = 1; i <= indexes.size(); i++) {
            index1 = indexes.get(indexes.size() - i);
            for (int j = 1; j < stString.length(); j++) {
                index0 = index1 - j;
                if (stString.charAt(index0) == ')') {
                    break;
                } else
                    index0--;
            }
            String newST = stString.substring(0, index0 + 1) + "I" + Integer.toString(i - 1) + stString.substring(index1);
            stString = newST;
        }
        System.out.println("AST String :" + stString);

        return stString;
    }


    // For GEMS: Check the RF_distance (trueGT, GEMS_GT) and (trueST, GEMS_GT)
    // Paras:  ifo.gtOutfitChecker("/Users/doriswang/Desktop/output/0001/10/50/", ifo.InputPath + "input/0001/", 10, 10);

    public void gtOutfitChecker(String iPath, String tPath, int lociNum, int fileNum, int seqLen) throws IOException {
        List<Tree> iGTs = loadIGTrees(iPath,lociNum, fileNum, seqLen);
        List<Tree> tGTs = new ArrayList<Tree>();
        String streeFile = tPath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        Tree trueST = Trees.readTree((String) stReader.readLine().trim());

        for (int ln = 0; ln < lociNum*fileNum; ln++) {
            tGTs.add(Trees.readTree(stReader.readLine().trim()));
        }
        stReader.close();
        double[] d = this.getDistances(iGTs,tGTs,fileNum*lociNum);

        BufferedWriter w = new BufferedWriter(new FileWriter(iPath  + "/distanceGEMS_GT.txt"));
        for(int i = 0;i<d.length;i++){
            w.write(i + ":" + d[i] + "\n");
            w.flush();
        }

        w.close();
        System.out.println("GEMS_Distance_GT:" + d[fileNum*lociNum]);

        double temp = 0.0;
        BufferedWriter w1 = new BufferedWriter(new FileWriter(iPath + "/" + lociNum + "/distanceGEMS_ST.txt"));
        for(int i = 0;i<lociNum*fileNum;i++){

            w1.write(i + ":" + this.getDistance(trueST,iGTs.get(i)) + "\n");
            temp += this.getDistance(trueST,iGTs.get(i));
            w1.flush();
        }
        w1.write(lociNum*fileNum + ":" + temp/(lociNum*fileNum));
        w1.flush();
        w1.close();

        System.out.println("GEMS_Distance_ST:" + temp/(lociNum*fileNum));


    }

    public double[] getDistances(List<Tree> l1, List<Tree> l2, int lociNum) throws IOException {
        Double dist = 0.0;
        double d[] = new double[lociNum + 1];

        for(int i = 0; i<l1.size(); i++){
            d[i] = this.getDistance(l1.get(i), l2.get(i));
            dist += this.getDistance(l1.get(i), l2.get(i));
        }
        d[lociNum] = dist/lociNum;
        return d;
    }


    // For GEMS: Check the RF_distance (trueGT, TF/R_GT) and (trueST, TF/R_GT)
    // Paras:  ifo.gtOutfitChecker("/Users/doriswang/Desktop/output/0001/10/50/", ifo.InputPath + "input/0001/", 10, 10);

    public void gtOverTFR(String iPath, String tPath, int lociNum, int fileNum, int seqLen) throws IOException {
        List<Tree> iGTs = new ArrayList<Tree>();

        List<Tree> tGTs = new ArrayList<Tree>();
        String streeFile = tPath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        Tree trueST = Trees.readTree((String) stReader.readLine().trim());

        for (int ln = 0; ln < lociNum*fileNum; ln++) {
            tGTs.add(Trees.readTree(stReader.readLine().trim()));
        }
        stReader.close();
        double[] d = this.getDistances(iGTs,tGTs,fileNum*lociNum);

        BufferedWriter w = new BufferedWriter(new FileWriter(iPath  + "/distanceTF/R_GT.txt"));
        for(int i = 0;i<d.length;i++){
            w.write(i + ":" + d[i] + "\n");
            w.flush();
        }

        w.close();
        System.out.println("TF/R_Distance_GT:" + d[fileNum*lociNum]);

        double temp = 0.0;
        BufferedWriter w1 = new BufferedWriter(new FileWriter(iPath + "/" + lociNum + "/distanceTF/R_ST.txt"));
        for(int i = 0;i<lociNum*fileNum;i++){

            w1.write(i + ":" + this.getDistance(trueST,iGTs.get(i)) + "\n");
            temp += this.getDistance(trueST,iGTs.get(i));
            w1.flush();
        }
        w1.write(lociNum*fileNum + ":" + temp/(lociNum*fileNum));
        w1.flush();
        w1.close();

        System.out.println("TF/R_Distance_ST:" + temp/(lociNum*fileNum));


    }

    // Compare with MCMC(Yang 2017)
    // Check the ratio: ( #true subtree appear / #inferred trees )
    public void MCMC_Compare() throws IOException {
        List iSTs = loadISTrees("/Users/doriswang/Desktop/a/0001/10/50/",2, 20);

        List subTrees = new ArrayList<Tree>();
//        subTrees.add(Trees.readTree("(1,2)"));
//        subTrees.add(Trees.readTree("(3,4)"));
//        subTrees.add(Trees.readTree("(5,6)"));
//        subTrees.add(Trees.readTree("(7,8)"));
//        subTrees.add(Trees.readTree("(9,10)"));
//        subTrees.add(Trees.readTree("(11,12)"));
//        subTrees.add(Trees.readTree("(13,14)"));
//        subTrees.add(Trees.readTree("(15,16)"));
//        subTrees.add(Trees.readTree("((1,2),(3,4))"));
//        subTrees.add(Trees.readTree("((5,6),(7,8))"));
//        subTrees.add(Trees.readTree("((9,10),(11,12))"));
//        subTrees.add(Trees.readTree("((13,14),(15,16))"));
//        subTrees.add(Trees.readTree("(((1,2),(3,4)),((5,6),(7,8)))"));
//        subTrees.add(Trees.readTree("(((9,10),(11,12)),((13,14),(15,16)))"));
//        subTrees.add(Trees.readTree("((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))"));
        //(1:1.0,(2:1.0,(3:1.0,(4:1.0,(5:1.0,(6:1.0,(7:1.0,(8:1.0,(9:1.0,(10:1.0,(11:1.0,(12:1.0,(13:1.0,(14:1.0,(15:1.0,16:1.0)I14:1.0)I13:1.0)I12:1.0)I11:1.0)I10:1.0)I9:1.0)I8:1.0)I7:1.0)I6:1.0)I5:1.0)I4:1.0)I3:1.0)I2:1.0)I1:1.0)I0;

        subTrees.add(Trees.readTree("(15,16)"));
        subTrees.add(Trees.readTree("((15,16),14)"));
        subTrees.add(Trees.readTree("(((15,16),14),13)"));
        subTrees.add(Trees.readTree("((((15,16),14),13),12)"));
        subTrees.add(Trees.readTree("(((((15,16),14),13),12),11)"));
        subTrees.add(Trees.readTree("((((((15,16),14),13),12),11),10)"));
        subTrees.add(Trees.readTree("(((((((15,16),14),13),12),11),10),9)"));
        subTrees.add(Trees.readTree("((((((((15,16),14),13),12),11),10),9),8)"));
        subTrees.add(Trees.readTree("(((((((((15,16),14),13),12),11),10),9),8),7)"));
        subTrees.add(Trees.readTree("((((((((((15,16),14),13),12),11),10),9),8),7),6)"));
        subTrees.add(Trees.readTree("(((((((((((15,16),14),13),12),11),10),9),8),7),6),5)"));
        subTrees.add(Trees.readTree("((((((((((((15,16),14),13),12),11),10),9),8),7),6),5),4)"));
        subTrees.add(Trees.readTree("((((((((((((15,16),14),13),12),11),10),9),8),7),6),5),3)"));
        subTrees.add(Trees.readTree("(((((((((((((15,16),14),13),12),11),10),9),8),7),6),5),3),2)"));
        subTrees.add(Trees.readTree("((((((((((((((15,16),14),13),12),11),10),9),8),7),6),5),3),2),1)"));


        double[] d = this.MCMCCompare(iSTs,subTrees);

        ////ifo.isTopoExist(t1,t2);
        // Syemmtric:
        //"((((16:1.0,15:1.0):1.0,(14:1.0,13:1.0):1.0):1.0,((12:1.0,11:1.0):1.0,(10:1.0,9:1.0):1.0):1.0):1.0,(((8:1.0,7:1.0):1.0,(6:1.0,5:1.0):1.0):1.0,((4:1.0,3:1.0):1.0,(2:1.0,1:1.0):1.0):1.0):1.0);"
        //Asymmetric:
        //
        Tree temp = Trees.readTree("(1:1.0,(2:1.0,(3:1.0,(4:1.0,(5:1.0,(6:1.0,(7:1.0,(8:1.0,(9:1.0,(10:1.0,(11:1.0,(12:1.0,(13:1.0,(14:1.0,(15:1.0,16:1.0)I14:1.0)I13:1.0)I12:1.0)I11:1.0)I10:1.0)I9:1.0)I8:1.0)I7:1.0)I6:1.0)I5:1.0)I4:1.0)I3:1.0)I2:1.0)I1:1.0)I0;");
        for(int i = 0;i<d.length;i++){
            System.out.println(i + " : " + d[i]);
        }
    }

    //load inferred species tree
    public static List<Tree> loadISTrees(String filePath,int lociNum, int treeNum) throws IOException {
        List<Tree> iST = new ArrayList<Tree>();
        for (int ln = 0; ln < treeNum; ln++) {
            String streeFile = "";
            if(lociNum == 2)
                 streeFile = filePath + ln + "/1000/50/GlobalBestTrees.txt";
            else if(lociNum == 10){
                streeFile = filePath + ln + "/1000/50/GlobalBestTrees.txt";
            }
            BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
            Tree ist = Trees.readTree(stReader.readLine().trim().split(": ")[1]);
            iST.add(ist);
        }
        return iST;
    }
    //load inferred gene tree
    public static List<Tree> loadIGTrees(String filePath,int lociNum, int treeNum, int seqLen) throws IOException {
        List<Tree> iGT = new ArrayList<Tree>();
        for (int ln = 0; ln < treeNum; ln++) {
            String streeFile = "";
            if(lociNum == 10 && seqLen==1000)
                streeFile = filePath + ln + "/" + seqLen + "/50/GlobalBestTrees.txt";
            else {
                streeFile = filePath + ln + "/" + seqLen + "/GlobalBestTrees.txt";

            }
            BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
            Tree igt = Trees.readTree(stReader.readLine().trim().split(": ")[1]);
            igt = Trees.readTree(stReader.readLine().trim().split(": ")[1]);
            iGT.add(igt);
            for(int i = 1;i<lociNum;i++){
                igt = Trees.readTree(stReader.readLine().trim());
                iGT.add(igt);
            }
            stReader.close();
        }
        return iGT;
    }

    public TNode getFatherNode(TNode n1, TNode n2){
        TNode big;
        TNode small;
        int l1 = n1.getLeafCount();
        int l2 = n2.getLeafCount();
        if(l1>=l2) {
            big = n1;
            small = n2;
        }
        else{
            big = n2;
            small = n1;
        }
        l1 = 0;
        boolean contain = true;
        List<String> smallList = new ArrayList<String>();
        List<String> bigList = new ArrayList<String>();
        Iterator it = small.getLeaves().iterator();
        while(it.hasNext()){
            String name = ((TNode)it.next()).getName();
            smallList.add(name);
        }
        it = big.getLeaves().iterator();
        while(it.hasNext()){
            String name = ((TNode)it.next()).getName();
            bigList.add(name);
        }
        for(int i = 0; i< smallList.size();i++){
            if(bigList.contains(smallList.get(i)))
                continue;
            else {
                contain = false;
                break;
            }
        }
        if(contain)
            return big;
        else
            return getFatherNode(big.getParent(),small);
    }

    public TNode getListFatherNode(List<TNode> nodes, Tree iST){

        TNode first = nodes.get(0);
        TNode tempFather = first;
        for(int i = 1;i<nodes.size();i++){
            TNode next = nodes.get(i);
            tempFather = getFatherNode(next, tempFather);
        }

        return tempFather;
    }


    public boolean isTopoExist(Tree iST, Tree subTree){

        String[] subNodes = subTree.getLeaves();
        List<TNode> nodes = new ArrayList<TNode>();
        for(int i =0;i<subNodes.length;i++){
            nodes.add(iST.getNode(subNodes[i]));
        }
        //find most smallest father node in iST for subNodes
        TNode t = getListFatherNode(nodes, iST);

        boolean isSameName = (subNodes.length==t.getLeafCount());
        boolean contain = true;
        if(!isSameName)
            return false;

        else {
            contain = true;
            List<String> list = new ArrayList<String>();
            Iterator it = t.getLeaves().iterator();
            while(it.hasNext()){
                String name = ((TNode)it.next()).getName();
                list.add(name);
            }

            for(int j = 0; j < list.size();j++){
                if(list.contains(subNodes[j]))
                    continue;
                else {
                    contain = false;
                    break;
                }
            }
        }
        return contain;
    }


    public boolean isSubTreeExist(TNode lNode, TNode rNode){

        TNode l = getNeighbor(lNode);
        TNode r = getNeighbor(rNode);

        if(l.isLeaf())
            return (l.getName()==r.getName());
        else if(r.isLeaf())
            return false;
        else
            return isSubTreeExist(l,r);
    }

    public TNode getNeighbor(TNode n){


        TNode parent = n.getParent();
        TNode temp;
        TNode neighbor = n;
        Iterator it = parent.getChildren().iterator();
        temp = (TNode)it.next();
        if(temp.getID() == n.getID())
            neighbor = (TNode)it.next();

        return neighbor;
    }

    public double[] MCMCCompare(List<Tree> iSTs, List<Tree> subTrees){
        double[] p = new double[subTrees.size()];
        int[][] counter = new int [iSTs.size()][subTrees.size()];
        for(int i = 0; i<iSTs.size();i++){
            Tree t = iSTs.get(i);
            for(int j = 0;j<subTrees.size();j++){
                if(isTopoExist(t,subTrees.get(j)))
                    counter[i][j] = 1;

                else
                    counter[i][j] = 0;
            }
        }
        double sum = 0;
        for(int j = 0; j < subTrees.size();j++){
            sum = 0;
            for(int i = 0;i<iSTs.size();i++){
                sum += counter[i][j];
            }

            p[j] = sum/iSTs.size();

        }

        return p;
    }
//        String streeFile = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST1/1/Tree/trueST.txt";
//        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
//        String trueST = (String) stReader.readLine().trim();
//        String seqPath = "/Users/doriswang/PhyloNet/Data/17-taxon/004/ST1/1/Seq/";
//        List trueGTS = new ArrayList<Tree>();
//        List trueSeq = new ArrayList<Alignment>();
//
//        //TODO : 6.5 : count bp Number to decide if continue
//        for (int ln = 0; ln < 30 ; ln++) {
//            String tree = stReader.readLine().trim();
//            trueGTS.add(tree);
//        }
//        stReader.close();
//        for (int ln = 0; ln < 30; ln++) {
//            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + ln + "_2000.nex"));
//            int j = 0;
//            while (j < 20) {
//                sReader.readLine();
//                j++;
//            }
//            String[] taxaName = new String[17];
//            String seqs = "";
//            String line = "";
//            Map<String, String> locus = new HashMap<String, String>();
//            for (int i = 0; i < 17; i++) {
//                line = sReader.readLine().trim();
//                String[] temp = line.split(" ");
//                taxaName[i] = temp[0];
//                if (temp[1].equals(""))
//                    seqs = temp[2];
//                else
//                    seqs = temp[1];
//                locus.put(taxaName[i], seqs.substring(0, 1000));
//            }
//            //locus.remove("O");
//            trueSeq.add(new Alignment(locus));
//
//            sReader.close();
//        }
//        ifo.initByRAxML(trueSeq,1000,17,30);
        //get16Aln(trueSeq, seqlens, 16,lociNum);


    public Map<String, Integer> loadSeqLens() throws IOException{
        Map<String, Integer> locusLen = new HashMap<String, Integer>();
        //<name, index>
        Map<String, String> locusName = new HashMap<String, String>();
        String inputFile = InputPath + "input/syr/syr011.txt";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedWriter bw = new BufferedWriter(new FileWriter( InputPath + "input/syr/locus_Mapping.txt"));
        BufferedWriter bwl = new BufferedWriter(new FileWriter(InputPath + "input/syr/locus_Len.txt"));

        for(int i = 0 ; i< 63; i++) {
            br.readLine();
        }
        for(int i = 0 ; i< 19; i++) {
            String temp = br.readLine().trim();
            String t[] = temp.split(" ");
            String lName = t[1];
            String[] tempLen = t[3].split("-");
            int start = Integer.valueOf(tempLen[0]);
            int end = Integer.valueOf((tempLen[1].split(";")[0]));
            System.out.print(String.valueOf(end - start+1)  + ",");
            locusLen.put(lName, end - start+1);
            locusName.put(lName, String.valueOf(i));
            bw.write(i + ":" + lName + "\n");
            bwl.write(i + ":" +  String.valueOf(end - start+1) + "\n");
            bw.flush();
        }
        bw.close();
        for(int i = 0 ; i< 14; i++) {
            br.readLine();
        }

        //Read Tree

        for(int k = 0 ; k < 19; k++) {
            String path = InputPath + "input/SYRData/" + k + "/";
            //ifo.isExitsPath(path);
            bw = new BufferedWriter(new FileWriter(path + "startTree.txt"));

            String temp = br.readLine().trim();
            int sIndex = temp.indexOf("]");
            String stString = temp.substring(sIndex+1);
            int index1 = 1;
            int index0 = 0;
            List<Integer> indexes = new ArrayList<Integer>();
            for (int i = 0; i < stString.length(); i++) {
                if (stString.charAt(i) == ')') {
                    indexes.add(i);
                }
            }
            for (int i = indexes.size() - 2; i >0; i--) {
                index0 = indexes.get(i);
                for (int j = 1; j < stString.length()-index0; j++) {
                    index1 = index0 + j;
                    if (stString.charAt(index1) == ':') {
                        break;
                    }
                }
                String newST = stString.substring(0, index0+1) + "I" + Integer.toString(i) + stString.substring(index1 );
                stString = newST;
            }
            bw.write(stString + "\n");
            bw.flush();
            bw.close();
        }
        return locusLen;
    }

    public Alignment loadLocus(int seqNum, int seqLength, int taxaNum, String inputFile) throws IOException, InterruptedException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        inputFile = inputFile + "seq_" + seqNum + "_" + seqLength + ".nex";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int j = 1;
        while (j < 21) {
            br.readLine();
            j++;
        }
        for (int i = 0; i < taxaNum; i++) {
            line = br.readLine().trim();
            String[] temp = line.split((" "));
            if (temp.length > 2) {
                temp[1] = temp[2];
            }
            String lociName = temp[0];
            String seq = "";
            seq = temp[1].substring(0, seqLength);
            locus.put(lociName, seq);
            //locus.remove("O");
        }
        Alignment aln = new Alignment(locus);
        isExitsPath(_RAxMLdir + seqNum );
        BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + seqNum + "/dna.phy"));
        phy.write(taxaNum + " " + seqLength + '\n');
        phy.flush();
        List<String> names = aln.getTaxaNames();
        Map<String, String> thisAln = aln.getAlignment();
        for (int i = 0; i < taxaNum; i++) {
            String name = names.get(i);
            String seq = thisAln.get(name);
            phy.write(name + '\t' + seq + '\n');
        }
        phy.flush();
        phy.close();

        return aln;
    }


    public Alignment loadRandomLocus(int seqNum, int seqLength, int taxaNum, String inputFile, int offset) throws IOException, InterruptedException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        inputFile = inputFile + "seq_" + String.valueOf(seqNum+offset) + "_" + seqLength + ".nex";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int j = 1;
        while (j < 21) {
            br.readLine();
            j++;
        }
        for (int i = 0; i < taxaNum; i++) {
            line = br.readLine().trim();
            String[] temp = line.split((" "));
            if (temp.length > 2) {
                temp[1] = temp[2];
            }
            String lociName = temp[0];
            String seq = "";
            seq = temp[1].substring(0, seqLength);
            locus.put(lociName, seq);
            //locus.remove("O");
        }
        Alignment aln = new Alignment(locus);
        isExitsPath(_RAxMLdir  + seqNum );
        BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir  + seqNum + "/dna.phy"));
        phy.write(taxaNum + " " + seqLength + '\n');
        phy.flush();
        List<String> names = aln.getTaxaNames();
        Map<String, String> thisAln = aln.getAlignment();
        for (int i = 0; i < taxaNum; i++) {
            String name = names.get(i);
            String seq = thisAln.get(name);
            phy.write(name + '\t' + seq + '\n');
        }
        phy.flush();
        phy.close();

        return aln;
    }

    //load seq data from raxml
    //TODO RAxML path
    public List<Alignment> loadSYRSeq48(int lociNum, String inputPath, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        List<Alignment> fullalns = new ArrayList<Alignment>();
        for(int l = 0; l<lociNum; l++){
            Map<String, String> locus = new HashMap<String, String>(); //20 taxa except 1
            Map<String, String> fullLocus = new HashMap<String, String>(); // 21 taxa include i(OG)
            String inputFile = inputPath + l + "/seq1.phy";
            String inputFile2 = inputPath + l + "/seq2.phy";
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            BufferedReader br2 = new BufferedReader(new FileReader(inputFile2));
            String line = "";
            String[] temp;
            String line2 = "";
            String[] temp2;
            line = br.readLine().trim();
            br2.readLine().trim();
            temp = line.split(" ");
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir  + l + "/dna.phy"));
            phy.write(taxaNum + " " + temp[1] + '\n');
            phy.flush();
            for (int i = 0; i < 26; i++) {
                line = br.readLine().trim();
                temp = line.split("\t");
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                //String seq = temp[1].substring(0, seqLength);
                locus.put(lociName, temp[1]);
                fullLocus.put(lociName, temp[1]);

                line2 = br2.readLine().trim();
                temp2 = line2.split("\t");
                if(temp2.length > 2)
                    temp2[1] = temp2[2];

                int tempIndex = Integer.valueOf(temp2[0].toString()) + 26;
                String lociName2 = String.valueOf(tempIndex);
                locus.put(lociName2,temp2[1]);
                fullLocus.put(lociName2, temp2[1]);
            }
            locus.remove("0");
            locus.remove("1");
            locus.remove("26");
            locus.remove("27");
            //locus.remove("21");
            br.close();
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullalns.add(fullAln);


            List<String> names = aln.getTaxaNames();
            Map<String, String> thisAln = aln.getAlignment();
            for (int i = 0; i < taxaNum; i++) {
                String name = names.get(i);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }

        return trueSeq;
    }

    //load seq data from raxml
    //TODO RAxML path
    public List<Alignment> loadSYRSeq(int lociNum, String inputPath, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
       List<Alignment> fullalns = new ArrayList<Alignment>();
        for(int l = 0; l<lociNum; l++){
            Map<String, String> locus = new HashMap<String, String>(); //20 taxa except 1
            Map<String, String> fullLocus = new HashMap<String, String>(); // 21 taxa include i(OG)
            String inputFile = inputPath + l + "/seq1.phy";
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            String line = "";
            String[] temp;
            line = br.readLine().trim();
            temp = line.split(" ");
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir  + l + "/dna.phy"));
            phy.write(taxaNum + " " + temp[1] + '\n');
            phy.flush();
            for (int i = 0; i < taxaNum+2; i++) {
                line = br.readLine().trim();
                temp = line.split(("\t"));
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                //String seq = temp[1].substring(0, seqLength);
                locus.put(lociName, temp[1]);
                fullLocus.put(lociName, temp[1]);
            }
            locus.remove("0");
            locus.remove("1");
            //locus.remove("3");
            //locus.remove("11");
            //locus.remove("21");
            br.close();
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullalns.add(fullAln);


            List<String> names = aln.getTaxaNames();
            Map<String, String> thisAln = aln.getAlignment();
            for (int i = 0; i < taxaNum; i++) {
                String name = names.get(i);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }

        return trueSeq;
    }

    //load seq data from raxml
    //TODO RAxML path
    public List<Alignment> loadSYRByRX(int lociNum, int[] seqLength, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>(); //20 taxa except 1
        Map<String, String> fullLocus = new HashMap<String, String>(); // 21 taxa include i(OG)
        List<Alignment> fullalns = new ArrayList<Alignment>();
        for(int l = 0; l<lociNum; l++){
            String inputFile = l + "/seq1.phy";
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            String line = "";
            String[] temp;
            line = br.readLine().trim();
            temp = line.split("\t");
            for (int i = 0; i < taxaNum+5; i++) {
                line = br.readLine().trim();
                temp = line.split(("\t"));
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                //String seq = temp[1].substring(0, seqLength);
                locus.put(lociName, temp[1]);
                fullLocus.put(lociName, temp[1]);
                locus.remove("0");
                locus.remove("1");
                locus.remove("3");
                locus.remove("11");
                locus.remove("21");
            }
            br.close();
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullalns.add(fullAln);
        }

        return trueSeq;
    }


    // set true alignment
    // return full alignment
    public List<Alignment> loadSYR(int indivNum, int[] seqLength, int taxaNum, List<Alignment> trueSeq) throws IOException {
        //"seq_" + tid + "_" + seqLen +".nex";scale + "/seq_" + i + "_" + len + ".nex"
        Map<String, String> locus = new HashMap<String, String>();
        Map<String, String> fullLocus = new HashMap<String, String>();
        List<Alignment> fullSeq = new ArrayList<Alignment>();
        String inputFile =InputPath + "input/syr/" + indivNum+ "/seq.txt";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String[] line = new String[seqLength.length];
        // List <taxon <seqs> >
        ArrayList<String[]> strArray = new ArrayList<String[]>();
        String[] taxaName = new String[taxaNum];
        for (int i = 0; i < taxaNum; i++) {
            line = new String[seqLength.length];
            String[] seqLines = br.readLine().trim().split(":");
            taxaName[i] = seqLines[0];
            String temp = seqLines[1];
            if (temp.length() < 8000)
                System.err.println("Loading failed at : " + i);
            int index = 0;
            for (int j = 0; j < seqLength.length; j++) {
                int len = seqLength[j];
                if (j == seqLength.length - 1)
                    line[j] = temp.substring(index);
                else {
                    line[j] = temp.substring(index, index + len);
                    index += len;
                }
            }
            strArray.add(line);
        }

        for(int i = 0; i < seqLength.length; i++){
            locus = new HashMap<String, String>();
            fullLocus = new HashMap<String, String>();
            for(int j = 0; j< taxaNum; j++){
                locus.put(taxaName[j], strArray.get(j)[i]);
                fullLocus.put(taxaName[j], strArray.get(j)[i]);
            }
            locus.remove("0");
            Alignment fullAln = new Alignment(fullLocus);
            Alignment aln = new Alignment(locus);
            trueSeq.add(aln);
            fullSeq.add(fullAln);
            BufferedWriter phy = new BufferedWriter(new FileWriter( _RAxMLdir + (i+19) + "/dna.phy" ));
            phy.write( taxaNum + " " + seqLength[i] + '\n');
            phy.flush();
            List<String> names = fullAln.getTaxaNames();
            Map<String, String> thisAln = fullAln.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        return fullSeq;
    }

    //Input: txt
    //Output: i/dna.phy


//  outgroup should be the child of root
    public String removeOutgroup(String tree) {
        //reroot
        // remove
        MutableTree thisTree = (MutableTree) Trees.readTree(tree.split(";")[0]);
        TMutableNode root = thisTree.getRoot();
        Iterator it = root.getChildren().iterator();
        TMutableNode otherChild = root;
        while (it.hasNext()) {
            TMutableNode child = (TMutableNode) it.next();
            if (child.getName().equals("O")) {
                if (it.hasNext())
                    otherChild = (TMutableNode) it.next();

                root.removeChild(child,false);
                thisTree.rerootTreeAtNode(otherChild);
                otherChild.removeChild(root,true);
                return thisTree.toString();
            }
            else
                otherChild = child;
        }
        return thisTree.toString();

    }

    //String seqFileName = "/Users/doriswang/PhyloNet/Data/IIG/seq/yeast_infer_network_from_multilocus_data_bayesian.nexus";
////t = 7 l = 106
    public List<Alignment> loadYData(String fileName, int taxaNum, int lociNum) throws IOException {
        List<Alignment> alns = new ArrayList<Alignment>();
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        int in = 0;
        //int[] lens = new int[106];
        while (in < 8) {
            br.readLine();
            in++;
        }
        int len = 0;
        for (int ln = 0; ln < lociNum; ln++) {
            Map<String, String> locus = new HashMap<String, String>();
            //List<Alignment> alns = new ArrayList<Alignment>();
            String line = "";
            String lociName = "";
            String seq = "";
            for (int i = 0; i < taxaNum; i++) {
                line = br.readLine().trim();

                String[] temp = line.split(" ");
                int index = 0;

                if (temp[index].length() > 1)
                    lociName = temp[0];
                else
                    lociName = temp[++index];
                if (temp[++index].length() > 1)
                    seq = temp[index];
                else
                    seq = temp[++index];
                locus.put(lociName, seq);
                //System.out.println("Locus_" + i +"'s length is " + seq.length());
                //len = seq.length();
            }
            //System.out.println("Locus_" + ln +"'s length is " + seq.length());
            System.out.print(" ," + seq.length());
            len += seq.length();
            //lens[lociNum] = seq.length();
            Alignment aln = new Alignment(locus);
            alns.add(aln);
            line = br.readLine();
            //if(line)
        }
        System.out.println("Ave length is " + len / 106);
        return alns;
    }


    //String seqFileName = "/Users/doriswang/PhyloNet/Data/IIG/seq/yeast_infer_network_from_multilocus_data_bayesian.nexus";
////t = 7 l = 106
    public List<Alignment> loadRData(String filePath, int repNum, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException {
        //#lociNum = 25
        // count bp Number to decide if continue
        for (int ln = 0; ln < lociNum; ln++) {
            String treeFile = filePath + "Rep" + repNum + ".gt" + ln + ".rose.tree.t";
            String seqFile = filePath + "Rep" + repNum + ".gt" + ln + ".rose.true.aln";
            BufferedReader tReader = new BufferedReader(new FileReader(treeFile));
            String tree = tReader.readLine().trim();
            trueGTS.add(tree);
            String[] lociName = new String[100];
            String[] seqs = new String[100];
            BufferedReader sReader = new BufferedReader(new FileReader(seqFile));
            String[] lines = sReader.readLine().trim().split(" ");
            //int tNum = Integer.getInteger(lines[0]);
            int tLength = Integer.valueOf(lines[1]);
            String line = "";
            for (int i = 0; i < 100; i++) {
                line = sReader.readLine().trim();
                String[] temp = line.split("    ");
                lociName[i] = temp[0];
                seqs[i] = temp[2].trim();
            }
            while (tLength != seqs[0].length()) {
                //SBF         GTGAATCAGGACGAAATTACTAGTAATACGCTAACAAGAATTCATGA--TATAGTGTCAA
                sReader.readLine();
                for (int i = 0; i < 100; i++) {
                    String l = sReader.readLine().trim();
                    seqs[i] += l;
                }
            }
            Map<String, String> locus = new HashMap<String, String>();
            for (int i = 0; i < 100; i++) {
                locus.put(lociName[i], seqs[i]);
            }
            trueSeq.add(new Alignment(locus));
        }
        //System.out.println("Ave length is " + len / 106);
        return trueSeq;
    }


//    //TODO: input: gtIndex: gtNum
//    public String loadRandomData (int stIndex, int seqlens, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException, ParseException, InterruptedException {
//        //#lociNum = 32
//        String streeFile = "/Users/doriswang/PhyloNet/Data/17-taxon/st/Rep" + stIndex +  "stree.txt";
//        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
//        String trueST = (String) stReader.readLine().trim();
//        String gtFile = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/Rep" + stIndex +  "gtrees.txt";
//        BufferedReader gtReader = new BufferedReader(new FileReader(streeFile));
//
//        for (int ln = 0; ln < lociNum ; ln++) {
//            String tree = stReader.readLine().trim();
//            trueGTS.add(tree);
//        }
//        stReader.close();
//
//        for (int ln = 0; ln < lociNum ; ln++) {
//            int tempLN = ln ;
//            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + tempLN + "_2000.nex"));
//            int j = 0;
//            while (j < 20) {
//                sReader.readLine();
//                j++;
//            }
//            String[] taxaName = new String[17];
//            String seqs = "";
//            String line = "";
//            Map<String, String> locus = new HashMap<String, String>();
//            for (int i = 0; i < 17; i++) {
//                line = sReader.readLine().trim();
//                String[] temp = line.split(" ");
//                taxaName[i] = temp[0];
//                if (temp[1].equals(""))
//                    seqs = temp[2];
//                else
//                    seqs = temp[1];
//                locus.put(taxaName[i], seqs.substring(0, seqlens));
//            }
//            //locus.remove("O");
//            trueSeq.add(new Alignment(locus));
//            sReader.close();
//        }
//        int taxaNum = 17;
//        int seqLen = 1000;
//        for(int i = 0; i<lociNum; i++){
//            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir + i + "/dna.phy"));
//            phy.write( taxaNum + " " + seqLen + '\n');
//            phy.flush();
//            Alignment a = trueSeq.get(i);
//            List<String> names = a.getTaxaNames();
//            Map<String, String> thisAln = a.getAlignment();
//            for (int j = 0; j < taxaNum; j++) {
//                String name = names.get(j);
//                String seq = thisAln.get(name);
//                phy.write(name + '\t' + seq + '\n');
//            }
//            phy.flush();
//            phy.close();
//        }
//        get16Aln(trueSeq, seqlens, 16,lociNum);
//        return trueST;
//    }

    //String seqFileName = "/Users/doriswang/PhyloNet/Data/17-taxon/32loci/Rep467gtrees";
////t = 16 + 1 l = 32 length = 2000
    //TODO 7.18 : only fit for 17 taxon now!! need to be extened
    public String load17Data (String filePath, int seqlens, int lociNum, List<Alignment> trueSeq, List<String> trueGTS) throws IOException, ParseException, InterruptedException {
        //#lociNum = 32
        String streeFile = filePath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        String trueST = (String) stReader.readLine().trim();
        String seqPath = filePath + "Seq/";
//        for(int t = 0;t<10;t++){
//            stReader.readLine().trim();
//        }
        // 6.5 : count bp Number to decide if continue
        stReader = new BufferedReader(new FileReader(filePath + "Tree/trueGTS.txt"));
        for (int ln = 0; ln < lociNum ; ln++) {
            String tree = stReader.readLine().trim();
            trueGTS.add(tree);
        }
        stReader.close();
        for (int ln = 0; ln < lociNum ; ln++) {
            int tempLN = ln ;
            //TODO different file name
            BufferedReader sReader = new BufferedReader(new FileReader(seqPath + "seq_" + tempLN + "_1000.nex"));
            int j = 0;
            while (j < 20) {
                sReader.readLine();
                j++;
            }
            String[] taxaName = new String[17];
            String seqs = "";
            String line = "";
            Map<String, String> locus = new HashMap<String, String>();
            for (int i = 0; i < 17; i++) {
                line = sReader.readLine().trim();
                String[] temp = line.split(" ");
                taxaName[i] = temp[0];
                if (temp[1].equals(""))
                    seqs = temp[2];
                else
                    seqs = temp[1];
                locus.put(taxaName[i], seqs.substring(0, seqlens));
            }
            //locus.remove("O");
            trueSeq.add(new Alignment(locus));
            sReader.close();
        }
        int taxaNum = 17;
        int seqLen = seqlens;
        for(int i = 0; i<lociNum; i++){
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir   + "/"+ i + "/dna.phy"));
            phy.write( taxaNum + " " + seqLen + '\n');
            phy.flush();
            Alignment a = trueSeq.get(i);
            List<String> names = a.getTaxaNames();
            Map<String, String> thisAln = a.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        get16Aln(trueSeq, seqlens, 16,lociNum);
        return trueST;
    }


    public  List<Tree> simulateGeneTrees(Network net,
                                               Map<String, List<String>> species2alleles, int numGT) {
        SimGTInNetworkByMS sim = new SimGTInNetworkByMS();
        return sim.generateGTs(Networks.readNetwork(net.toString()), species2alleles, numGT, _msdir);
    }

    //seq-gen file
    public List<Alignment> loadAlnSQ(String inputFile, int seqLength, int taxaNum, int loci) throws IOException {
        Map<String, String> locus = new HashMap<String, String>();
        List<Alignment> alns = new ArrayList<Alignment>();
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = "";
        int lineNum = 1;
        while (lineNum < 20) {
            br.readLine();
            lineNum++;
        }
        for (int j = 0; j < loci; j++) {
            if (j > 0) {
                int inter = 1;
                while (inter < 8) {
                    br.readLine();
                    inter++;
                }
            }
            for (int i = 0; i < taxaNum; i++) {
                line = br.readLine().trim();

                String[] temp = line.split(" ");
                if (temp.length > 2) {
                    temp[1] = temp[2];
                }
                String lociName = temp[0];
                String seq = "";
                seq = temp[1].substring(0, seqLength);
                locus.put(lociName, seq);
            }
            Alignment aln = new Alignment(locus);
            alns.add(aln);
        }

        return alns;
    }


    public static ArrayList<File> getListFiles(Object obj) {
        File directory = null;
        if (obj instanceof File) {
            directory = (File) obj;
        } else {
            directory = new File(obj.toString());
        }
        ArrayList<File> files = new ArrayList<File>();
        if (directory.isFile()) {
            files.add(directory);
            return files;
        } else if (directory.isDirectory()) {
            File[] fileArr = directory.listFiles();
            for (int i = 0; i < fileArr.length; i++) {
                File fileOne = fileArr[i];
                files.addAll(getListFiles(fileOne));
            }
        }
        return files;
    }

    public List<UltrametricTree> simGTByUPGMA(Alignment seq) {
        //BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));

        UPGMATree upgma = new UPGMATree(new JCDistance(seq.getAlignment()));
        UltrametricTree template = upgma.getUltrametricTree();
        //System.out.println("Init gene tree： " + template.getTree().toNewick());
        List<UltrametricTree> geneTrees = new ArrayList<>();
        geneTrees.add(upgma.getUltrametricTree());
        return geneTrees;
    }

    public List<UltrametricTree> simGTSByUPGMA(List<Alignment> trueAlns) throws IOException {
        List<UltrametricTree> geneTrees = new ArrayList<UltrametricTree>();
        for (int i = 0; i < trueAlns.size(); i++) {
            geneTrees.addAll(simGTByUPGMA((Alignment) trueAlns.get(i)));
        }
        return geneTrees;
    }


    protected Tree getGLASSTree(List<UltrametricTree> gtList, Map<String, List<String>> species2alleles) {
        GLASSInference glassInference = new GLASSInference();
        Tree inferredST = glassInference.inferSpeciesTree(getTreeByUTree(gtList));
        return inferredST;
    }

    //
    public double getAverage(double[] array){
        int sum = 0;
        for(int i = 0;i < array.length;i++){
            sum += array[i];
        }
        return (double)(sum / array.length);
    }

    // P(gt|ST)
    //if not return -> add st height on all leaves
    public List<Double> getGTSLLBySTYF(List<Tree> geneTrees, Network<NetNodeInfo> currentST) {
        List<Double> gProST = new ArrayList<>();
        double[] result = new double[geneTrees.size()];
        Network<NetNodeInfo> ultraST = getUltraST(currentST.toString());
        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(ultraST, geneTrees, null);
        g.calculateGTDistribution(result);
        for (int i = 0; i < result.length; i++)
            gProST.add(Math.log(result[i]));

        return gProST;
    }

//    // input: gtLL( P(seq|topos), topos (refined by RAxML), st
//    //output:
//    //if not return -> add st height on all leaves
//    public List<Double> getRGTSLLBySTYF(List<Double> gtLL, List<Tree> topos, Network<NetNodeInfo> currentST) {
//        List<Double> gProST = new ArrayList<>();
//        double[] result = new double[geneTrees.size()];
//        Network<NetNodeInfo> ultraST = getUltraST(currentST.toString());
//        GeneTreeWithBranchLengthProbabilityYF g = new GeneTreeWithBranchLengthProbabilityYF(ultraST, geneTrees, null);
//        g.calculateGTDistribution(result);
//        for (int i = 0; i < result.length; i++)
//            gProST.add(Math.log(result[i]));
//
//        return gProST;
//    }

    //
    public double getStandardDevition(double[] array){
        double sum = 0;
        double ave = getAverage(array);
        for(int i = 0;i < array.length;i++){
            sum += Math.sqrt(((double)array[i] - ave) * (array[i] - ave));
        }
        return (sum / (array.length - 1));
    }
    public BniNetwork<NetNodeInfo> inferSTByGLASS(List<UltrametricTree> geneTrees, double halfTheta) throws IOException, ParseException {
        Tree speciesTree = getGLASSTree(geneTrees, null);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils._POP_SIZE_MEAN = 0.036;
        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(geneTrees, null, null);
        UltrametricNetwork stWithPara = new UltrametricNetwork(speciesTree.toNewick(), geneTrees, null);
        //stWithPara.initNetHeights(4*scale,constraints);
        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
        BniNetwork<NetNodeInfo> st = (BniNetwork) stWithPara.getNetwork();
        Networks.autoLabelNodes(st);
        return st;
    }

    public Tree inferSTByGLASS(List<Tree> geneTrees) throws IOException, ParseException {
        GLASSInference glassInference = new GLASSInference();
        Tree inferredST = glassInference.inferSpeciesTree(geneTrees);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils._POP_SIZE_MEAN = 0.001;
        List<UltrametricTree> uList = getUTreeByTree(geneTrees);
        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(uList, null, null);
        UltrametricNetwork stWithPara = new UltrametricNetwork(inferredST.toNewick(), uList, null);
        //stWithPara.initNetHeights(4*scale,constraints);
        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
        BniNetwork<NetNodeInfo> st = (BniNetwork) stWithPara.getNetwork();
        Networks.autoLabelNodes(st);
        Tree glassTree = Trees.readTree(st.toString());
        return glassTree;
    }

    //infer ST by GTS by ASTRAL
    public BniNetwork<NetNodeInfo> inferSTByAST(List<UltrametricTree> geneTrees, double popSize) throws IOException, ParseException {
        String command = "java -jar " + _ASTRALdir + "astral.4.11.1.jar -i " + _ASTRALdir + "test_data/tempIn.tre -o " + _ASTRALdir + "test_data/testOut.tre"; ///Users/doriswang/PhyloNet/tools/Astral/test_data/tempIn.tre -o /Users/doriswang/PhyloNet/tools/Astral/test_data/testOut.tre";

        String cmdFile = _ASTRALdir + "tempCMD.sh";
        BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
        cmd.write(command + '\n');
        cmd.flush();
        cmd.close();
        String inputFile = _ASTRALdir + "test_data/tempIn.tre";
        BufferedWriter in = new BufferedWriter(new FileWriter(inputFile));
        for (int i = 0; i < geneTrees.size(); i++) {
            in.write(geneTrees.get(i).toString() + "\n");
        }
        in.flush();
        in.close();
        ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
        pb.redirectErrorStream(true);
        try {
            Process proc = pb.start();
            try {
                proc.waitFor();
            } catch (InterruptedException ex) {
                System.out.println(ex.getMessage());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        String outputFile = _ASTRALdir + "test_data/testOut.tre";
        BufferedReader out = new BufferedReader(new FileReader(outputFile));
        String stString = out.readLine().trim();
        BniNetwork<NetNodeInfo> st1 = new BniNetwork<>();
        int index0 = 0;
        int index1 = 1;
        int count = 0;
        List<Integer> indexes = new ArrayList<Integer>();
        for (int i = 0; i < stString.length(); i++) {
            if (stString.charAt(index1) == ':') {
                indexes.add(i);
            }
        }
        for (int i = 1; i <= indexes.size(); i++) {
            index1 = indexes.get(indexes.size() - i);
            for (int j = 1; j < stString.length(); j++) {
                index0 = index1 - j;
                if (stString.charAt(index0) == ')') {
                    break;
                }
            }
            String newST = stString.substring(0, index0) + "I" + Integer.toString(i - 1) + stString.substring(index1 + 1);
            stString = newST;
        }
//        ArrayList<Double> paras = new ArrayList<Double>();
//        String topo = "";
//        String firstNode = "";
//        //String pattern1 = "^\\d+";
//        String pattern = "\\d+(\\.\\d+)?";
//        Pattern p = Pattern.compile(pattern);
//        String[] s1 = stString.split(":");
//        for(int i = 0; i<s1.length; i++){
//            Matcher m = p.matcher(s1[i]);
//            boolean hasStart = false;
//            int firstEnd = -1;
//            String temp = "";
//            while(m.find()){
//                int start = m.start();
//                int end = m.end();
//                if(start==0) {
//                    temp += s1[i].substring(end);
//                    hasStart = true;
//                    firstEnd = end;
//                    continue;
//                }
//                if(end==s1[i].length()) {
//                    if(hasStart)
//                        temp = s1[i].substring(firstEnd,start);
//                    else
//                        temp = s1[i].substring(0, start);
//                }
//            }
//            topo+=temp;
//        }
//
//        BniNetwork st = (BniNetwork) Networks.readNetwork(topo);
//        Networks.autoLabelNodes(st);
//        Utils._ESTIMATE_POP_SIZE = false;
//        Utils._CONST_POP_SIZE = true;
//        Utils._POP_SIZE_MEAN = popSize;
//        Map<String, Double> constraints = TemporalConstraints.getTemporalConstraints(geneTrees, null, null);
//        UltrametricNetwork stWithPara = new UltrametricNetwork(st.toString(), geneTrees, null);
//        //stWithPara.initNetHeights(4*scale,constraints);
//        //Networks.autoLabelNodes((Network)stWithPara.getNetwork());
//        BniNetwork<NetNodeInfo> st1 = (BniNetwork) stWithPara.getNetwork();
//        Networks.autoLabelNodes(st1);
//        return st1;

        return st1;
    }

    //./raxmlHPC -M -m GTRGAMMA -p 12345 -q simpleDNApartition.txt -s dna.phy -n T22
    public List<Tree> initByRAxML(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            System.out.println("RAXML No" + i + "Start...");
            isExitsPath(_RAxMLdir + i);
            String rm = "rm *.T" + i + '\n';
            //./raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 10 -s 0/dna.phy -n T0
            String raxml = _RAxMLdir + "raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 30 " + "-s " + i + "/dna.phy -n T" + i;
            String cmdFile = _RAxMLdir + i + "/tempCMD.sh";
            BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
            String tempR = _RAxMLdir ;
            cmd.write("cd " + tempR.substring(0, tempR.length() - 1) + '\n');
            cmd.write(rm);
            cmd.flush();
            cmd.write(raxml);
            cmd.flush();
            cmd.close();

            ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
            pb.redirectErrorStream(true);
            try {
                Process proc = pb.start();
                try {
                    proc.waitFor();
                } catch (InterruptedException ex) {
                    System.out.println(ex.getMessage());
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < lociNum; i++) {
            String outputFile = _RAxMLdir + "RAxML_bestTree.T" + i;
            BufferedReader gtReader = new BufferedReader(new FileReader(outputFile));
            String gtString = gtReader.readLine().trim();
            gts.add(Trees.readTree(gtString));
            gtReader.close();
        }
        return gts;
    }
//FastTree -gtr -nt < alignment.file > tree_file
    public List<Tree> initFasttree(List<Alignment> aln, int[] seqLen, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in

        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            //System.out.println("No" + i + "Start...");
            isExitsPath(_fastTree + i );
            BufferedWriter bw = new BufferedWriter(new FileWriter(_fastTree + "/" + i + ".phy"));
            Alignment temp = aln.get(i);
            Map<String, String> thisAln = temp.getAlignment();
            int taxaNum = temp.getTaxonSize();
            bw.write(taxaNum + " " + seqLen[i] + '\n');
            for (int j = 0; j < taxaNum; j++) {
                String name = temp.getTaxaNames().get(j);
                String seq = thisAln.get(name);
                bw.write(name + ' ' + seq + '\n');
            }
            bw.flush();
            bw.close();
            isExitsPath(_fastTree );
            String rm = "rm *.T" + i + '\n';
            //./raxmlHPC-PTHREADS -M -m GTRGAMMA -p 12345 -# 10 -s 0/dna.phy -n T0
            String command = _fastTree + "FastTree -gtr -nt < " + i + ".phy > " + i + ".T";
            String cmdFile = _fastTree + "/FastTree" + i + ".sh";
            BufferedWriter cmd = new BufferedWriter(new FileWriter(cmdFile));
            cmd.write("cd " + _fastTree.substring(0, _fastTree.length() - 1) + '\n');
            cmd.write(rm);
            cmd.flush();
            cmd.write(command);
            cmd.flush();
            cmd.close();

            ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
            pb.redirectErrorStream(true);
            try {
                Process proc = pb.start();
                try {
                    proc.waitFor();
                } catch (InterruptedException ex) {
                    System.out.println(ex.getMessage());
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        for (int i = 0; i < lociNum; i++) {
            String outputFile = _fastTree + i + ".T";
            BufferedReader out = new BufferedReader(new FileReader(outputFile));
            //String gtString = gtReader.readLine().trim();
            String stString = out.readLine().trim();
            int index0 = 0;
            int index1 = 1;
            for(int j = 1; j< stString.length()-1;j++){
                if(stString.charAt(j)!=')'){
                    continue;
                }
                else if(stString.charAt(j+1)==':'){
                    continue;
                }
                else{
                    int offset = 1;
                    boolean find = false;
                    for(offset=1;offset<stString.length()-j;offset++){
                        if(stString.charAt(offset+j)==':') {

                            find = true;
                            break;
                        }
                        else
                            continue;
                    }
                    if(find){
                        stString = stString.substring(0,j+1) + stString.substring(offset+j);
                    }
                }
            }
            gts.add(Trees.readTree(stString));
            //System.out.println(stString);
            out.close();
        }
        return gts;
    }
    //TODO
    public String[] getRAxMLStartT(int lociNum) throws IOException, ParseException, InterruptedException {
        //List<Tree> gts = new ArrayList<Tree>();
        String[] gts = new String[lociNum];
        for (int i = 0; i < lociNum; i++) {
            String tree = _RAxMLdir + "/RAxML_Init/RAxML_bestTree.T" + i;
            String llFile = _RAxMLdir + "/RAxML_Init/RAxML_info.T" + i;
            BufferedReader gtReader = new BufferedReader(new FileReader(tree));
            String gtString = gtReader.readLine().trim();
            gtReader.close();
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            double ll = 0.0;
            String llString = "";
            while(llString!=null){
                llString = llReader.readLine();
                if(llString.length()<20||llString.isEmpty())
                    continue;
                if(llString.substring(0,5).equals("Final")){
                    String[] result = llString.split(" ");
                    ll = Double.valueOf(result[6]);
                    break;
                }
            }
            gts[i] = gtString + " " + ll;
            llReader.close();
        }
        return gts;
    }


        //operator.getUpdateSh(lociNum, msSize);
    public void getUpdateSh(int lociNum, int msSize) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            String seqPath;
            for(int j = 0; j<msSize; j++){
                seqPath = _RAxMLdir  + i + "/";
                String tName = "L" + String.valueOf(i) + "T" + String.valueOf(j); // L_i T_j
                BufferedWriter phy = new BufferedWriter(new FileWriter(seqPath + "updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh"));
                System.out.println(seqPath + "updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh");
                phy.write("cd " + _RAxMLdir + '\n');
                phy.flush();
                phy.write("rm *." + tName + '\n');
                phy.flush();
                phy.write(_RAxMLdir    + "raxmlHPC-PTHREADS  -f e -t topo.T" + String.valueOf(j) + " -m GTRGAMMA -s " + String.valueOf(i) + "/dna.phy -n " + tName + '\n');
                phy.flush();
                phy.close();
            }

        }
        //return gts;
    }


    //operator.getUpdateSh(lociNum, msSize);
    public void changeUpdateSh(int lociNum, int msSize) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {

            String seqPath = _RAxMLdir + i + "/";
            for(int j = 0; j<msSize; j++){
                String tName = "L" + String.valueOf(i) + "T" + String.valueOf(j); // L_i T_j
                BufferedWriter phy = new BufferedWriter(new FileWriter(seqPath + "updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh"));
                //phy.write("cd " + "/home/yw58/IIG/tools/RAxML/1/" + '\n');
                phy.write("cd " + _RAxMLdir + '\n');

                phy.flush();
                phy.write("rm *." + tName + '\n');
                phy.flush();

                phy.write(_RAxMLdir + "raxmlHPC-PTHREADS  -f e -t topo.T" + String.valueOf(j) + " -m GTRGAMMA -s " + String.valueOf(i) + "/dna.phy -n " + tName + '\n');
                //phy.write("/home/yw58/IIG/tools/RAxML/1/" + "raxmlHPC-PTHREADS  -f e -t topo.T" + String.valueOf(j) + " -m GTRGAMMA -s " + String.valueOf(i) + "/dna.phy -n " + tName + '\n');
                phy.flush();
                phy.close();
            }

        }
        //return gts;
    }

    //write alignment with 16 taxa for RAxML
    public List<Tree> get16Aln(List<Alignment> aln, int seqLen, int taxaNum, int lociNum) throws IOException, ParseException, InterruptedException {
        //writeSeq + partition     rm + raxml    read in
        List<Tree> gts = new ArrayList<Tree>();
        for (int i = 0; i < lociNum; i++) {
            isExitsPath(_RAxMLdir +i);
            BufferedWriter phy = new BufferedWriter(new FileWriter(_RAxMLdir  + i + "/dnaNoOG.phy"));
            phy.write(taxaNum + " " + seqLen + '\n');
            phy.flush();
            Alignment a = aln.get(i);
            List<String> names = a.getTaxaNames();
            Map<String, String> thisAln = a.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                if (name == "O")
                    continue;
                String seq = thisAln.get(name);
                phy.write(name + '\t' + seq + '\n');
            }
            phy.flush();
            phy.close();
        }
        return gts;
    }

    //Undo: false -> *halftheta    true -> /halfTheta
    public static Network<NetNodeInfo> scaleSTNet(Network<NetNodeInfo> currentST, double ratio, boolean undo) {
        //Tree newT = new STITree<>();
        Iterator itNode = currentST.getTreeNodes().iterator();
        if (undo == false) {
            while (itNode.hasNext()) {
                NetNode n = (NetNode) itNode.next();
                if (n == (NetNode) currentST.getRoot())
                    continue;
                NetNode p = (NetNode) n.getParents().iterator().next();
                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(p, n.getParentDistance(p) * ratio); // height = height * ratio
                }
            }
        } else {
            while (itNode.hasNext()) {
                NetNode n = (NetNode) itNode.next();
                if (n == (NetNode) currentST.getRoot())
                    continue;
                NetNode p = (NetNode) n.getParents().iterator().next();
                if (n.getParentDistance(p) != TMutableNode.NO_DISTANCE) {
                    n.setParentDistance(p, n.getParentDistance(p) / ratio); // height = height / ratio
                }
            }
        }
        return currentST;
    }

    public boolean hasCheckedBL(BniNetNode node, Network st) {
        boolean hasChecked = true;
        while (!node.isLeaf()) {
            BniNetNode child = (BniNetNode) node.getChildren().iterator().next();
            if (child.getParentDistance(node) == Double.NEGATIVE_INFINITY)
                hasChecked = false;
            node = child;
        }
        return hasChecked;
    }

    //get distance(a,b) where common parent is cRoot on one gt
    public double getLocalDist(TNode left, TNode right, TNode cRoot, Tree gt) {
        double lDistance = left.getParentDistance();
        double rDistance = right.getParentDistance();

        TNode lParent = left.getParent();
        TNode rParent = right.getParent();
        while (lParent != cRoot) {
            lDistance += lParent.getParentDistance();
            lParent = lParent.getParent();
        }
        while (rParent != cRoot) {
            rDistance += rParent.getParentDistance();
            rParent = rParent.getParent();
        }
        return lDistance + rDistance;
    }
//    //p is child or child's ancestor
//    public boolean isAncestor(TNode parent, TNode child, Stack<TNode> s){
//        if(parent.equals(child))
//            return true;
//        else{
//            if(parent.isLeaf())
//                return false;
//            else {
//                Iterator it = parent.getChildren().iterator();
//                TNode temp1 = (TNode)it.next();
//                TNode temp2 = (TNode)it.next();
//                if(temp1.equals(child)||temp2.equals(child))
//                    return true;
//                else
//                    return isAncestor(temp1,child,s)|| isAncestor(temp2,child,s);
//            }
//        }
//    }


    //Input: true gene trees ; ms sampled trees for iteration
    //Output: the best distance IIG may get for each locus in each iteration
    public void getMSDist(List<Tree> trueTopos, List<Tree> msTrees, double[][] msTreeDist, int itNum){
        for(int i = 0; i<trueTopos.size();i++){
            Tree thistree = trueTopos.get(i);
            double minD = 1;
            for(int j = 0; j<msTrees.size();j++){
                Tree msT = msTrees.get(j);
                double d = getDistance(thistree,msT);
                if(d<minD)
                    minD = d;
            }
            msTreeDist[itNum][i] = minD;
        }
    }


    //get Lowest Common Parent on gt for computing distance(a,b)
    public TNode getLCP(TNode root, TNode p, TNode q) {
        //System.out.println("getLCP" + p.getName() + " -- " + q.getName());
        if (root == null)
            return null;
        if (root == p || root == q)
            return root;
        if (root.isLeaf())
            return null;
        else {
            Iterator it = root.getChildren().iterator();
            TNode l = getLCP((TNode) it.next(), p, q);
            TNode r = getLCP((TNode) it.next(), p, q);
            if (l != null && r != null) {
                return root;
            } else if (l == null && r == null) {
                return null;
            } else {
                return l == null ? r : l;
            }
        }
    }

    //getInternalHeight(sib, st, gts);
    public double getInternalHeight(NetNode thisNode, Network st, Map<String, Double> heightList) {
        double height = Double.NEGATIVE_INFINITY;
        Iterator children = thisNode.getChildren().iterator();
        while (children.hasNext()) {
            NetNode c = (NetNode) children.next();
            Double cHeight = 0.0;
            if(!heightList.containsKey(c.getName()))
                cHeight = getInternalHeight(c, st, heightList);
            if (height < cHeight + c.getParentDistance(thisNode)) {
                    height = heightList.get(c.getName()) + c.getParentDistance(thisNode);
                    heightList.put(thisNode.getName(), height);
                }
        }
        return height;

    }

    // return Map<#level, nodesList>
    public Map<Integer, List<String>> getLevels(Network st) {
        Map levelMap = new HashMap<Integer, List<String>>();
        LinkedList q = new LinkedList<NetNode>();
        q.offer(st.getRoot());
        int cur, last;
        int level = 1;
        NetNode current = new BniNetNode();
        // Map<Integer, String> map = new HashMap<Integer, String>();
        while (!q.isEmpty()) {
            cur = 0;
            last = q.size();
            while (cur < last) {
                current = (NetNode) q.poll();
                cur++;
                if (current.isLeaf()) {
                    if (levelMap.containsKey(level)) {
                        List nodes = (List) levelMap.get(level);
                        nodes.add(current.getName());
                    } else {
                        List<String> nodes = new ArrayList<String>();
                        nodes.add(current.getName());
                        levelMap.put(level, nodes);
                    }
                } else {
                    Iterator children = current.getChildren().iterator();
                    while (children.hasNext())
                        q.offer(children.next());
                }
            }
            level++;
        }
        return levelMap;
    }

    public Network<NetNodeInfo> getUltraST(String stString) {
        Network<NetNodeInfo> st = Networks.readNetwork(stString);
        List<NetNode> visit = new ArrayList<NetNode>();
        Map<NetNode,Double> heights = new HashMap<NetNode,Double>();
        Iterator leafIt = st.getLeaves().iterator();
        Double netHeight = 0.0;
        while(leafIt.hasNext()){
            NetNode temp = (NetNode)leafIt.next();
            if(visit.contains(temp))
                continue;
            else{
                visit.add(temp);
                double tempHeight = getLeafHeight(temp,st);
                if(tempHeight>netHeight)
                    netHeight = tempHeight;
                heights.put(temp,tempHeight);
                NetNode p = (NetNode) temp.getParents().iterator().next();
                NetNode sib = (NetNode) p.getChildren().iterator().next();
                if(sib.equals(temp))
                    sib = (NetNode) p.getChildren().iterator().next();
                if(sib.isLeaf()){
                    visit.add(sib);
                    heights.put(sib,tempHeight);
                }

            }
        }
        Iterator leaves = heights.keySet().iterator();
        while(leaves.hasNext()){
            NetNode thisNode = (NetNode)leaves.next();
            NetNode parent = (NetNode)thisNode.getParents().iterator().next();
            Double currentH = heights.get(thisNode);
            Double dist = thisNode.getParentDistance(parent);
            thisNode.setParentDistance(parent,dist + netHeight-currentH);

        }
        return st;
    }

    //    public double getLowestHeight(NetNode n, Network st){
//        double h = 0.0;
//        if(n.isLeaf())
//            return h;
//        Stack s = new Stack();
//        Iterator tempIt = n.getChildren().iterator();
//        s.add(n);
//        Map<String,Double> heights = new HashMap<String,Double>();
//        while(!s.isEmpty()){
//            NetNode thisNode = (NetNode) s.pop();
//            if(thisNode.isLeaf()){
//                heights.add
//            }
//        }
//        return h;
//    }
    public boolean isExitsPath(String path) throws InterruptedException {
        File file = new File(path.toString());
        if (!file.exists()) {
            file.mkdirs();
            System.out.println("build dir：" + path.toString());
            //Thread.sleep(1500);
        }
        File file1 = new File(path.toString());
        if (!file1.exists()) {
            return true;
        } else {
            return false;
        }
    }

    public void updateMSTopos(List<Tree> topos) throws IOException {

            String seqPath = _RAxMLdir ;
            for(int j = 0; j<topos.size(); j++){
                BufferedWriter topow = new BufferedWriter(new FileWriter(seqPath + "topo.T" + String.valueOf(j)));
                topow.write(topos.get(j).toString());
                topow.flush();
                topow.close();
            }

    }


    public void runUpdateShell(int msSize, int lociNum) throws IOException {

        //String seqPath = _RAxMLdir;
        for (int i = 0; i < lociNum; i++) {
            String seqPath = _RAxMLdir  + i + "/";
            for(int j = 0; j<msSize; j++){
                //updateL0T1.sh
                String cmdFile = _RAxMLdir  + i + "/updateL" + String.valueOf(i) + "T" + String.valueOf(j) + ".sh";
                ProcessBuilder pb = new ProcessBuilder("/bin/bash", cmdFile);
                pb.redirectErrorStream(true);
                try {
                    Process proc = pb.start();
                    try {
                        proc.waitFor();
                    } catch (InterruptedException ex) {
                        System.out.println(ex.getMessage());
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }


    //Input: #locus, refineSize
    //Output: double[] ll List for Locus_i
    public double[] getLocusLL(int locusNum, int msSize) throws IOException {
        double[] ll = new double[msSize];
        //String outputFile = _RAxMLdir + "RAxML_result.TEST" + i;
        for(int j = 0; j<msSize; j++){
            String llFile = _RAxMLdir  + "/RAxML_log.L" + String.valueOf(locusNum) + "T" + String.valueOf(j);
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            String[] llString = llReader.readLine().trim().split(" ");
            double tempLL = Double.valueOf(llString[1]);
            ll[j] = tempLL;
            llReader.close();
        }
        return ll;
    }


    //Input: trees[], ll[] ,lociNum
    //Output: double[] ll, trees[]
    public double[] getInitTree(String[] trees, double[] ll, int lociNum) throws IOException {

        for(int j = 0; j<lociNum; j++){
            String tFile = _RAxMLdir  + "RAxML_bestTree.T" + String.valueOf(j);
            BufferedReader tReader = new BufferedReader(new FileReader(tFile));
            trees[j] = tReader.readLine().trim();
            tReader.close();
            String llFile = _RAxMLdir + "RAxML_info.T" + String.valueOf(j);
            BufferedReader llReader = new BufferedReader(new FileReader(llFile));
            while(true){
                String temp = llReader.readLine().trim();
                if(temp.startsWith("Final")){
                    String temps[] = temp.split(" ");
                    temp = temps[temps.length-1];
                    ll[j] = Double.valueOf(temp);
                    break;
                }
            }
            llReader.close();
            //System.out.println(trees[j] + "  " + ll[j]);
        }
        return ll;
    }

    //Input: #locus, #bestMSNum
    //Output: RAxML tree
    public Tree getbestGTi(int locusNum, int msTreeNum) throws IOException {
        String gtFile = _RAxMLdir  + "RAxML_result.L" + String.valueOf(locusNum) + "T" + String.valueOf(msTreeNum);
        BufferedReader gtReader = new BufferedReader(new FileReader(gtFile));
        String gtString = gtReader.readLine().trim();
        Tree bestGT = Trees.readTree(gtString);
        gtReader.close();
        return bestGT;
    }


    //TODO check
    public static List<UltrametricTree> scaleGTS(List<UltrametricTree> uList, double ratio, boolean undo) {
        Iterator gList = uList.iterator();
        while (gList.hasNext()) {
            UltrametricTree g = (UltrametricTree) gList.next();
            Iterator itNode = g.getTree().getNodes().iterator();
            if (undo == false) {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) g.getTree().getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() * ratio); // height = height * ratio
                    }
                }
            } else {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) g.getTree().getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() / ratio); // height = height / ratio
                    }
                }
            }

        }
        return uList;
    }

    public boolean isRooted(Tree t) {
        boolean flag = true;
        TNode root = t.getRoot();
        if(root.getChildCount()>2)
            flag = false;
        Iterator cIt = root.getChildren().iterator();
        TNode c1 = (TNode) cIt.next();
        TNode c2 = (TNode) cIt.next();
        if(c1.getName().equals("O")||c2.getName().equals("O"))
            flag =  true;
        else
            flag = false;
        return flag;
    }

    //scale mtDNA
    public static List<Tree> scaleGTList(boolean multiply, List<Tree> tList, double ratio) {
        List<Tree> newGTS = new ArrayList<Tree>();
        for (int i = 0; i < tList.size(); i++) {
            newGTS.add(scaleGT(tList.get(i), ratio, multiply));
        }
        return tList;
    }

    public boolean equalTree(Tree t1, Tree t2){
        TNode n1 = t1.getRoot();
        TNode n2 = t2.getRoot();
        return ifSameSubTree(n1,n2);
    }


    public boolean ifSameLeaf(TNode n1 , TNode n2){
        return n1.getName().equals(n2.getName());
    }

    public boolean ifSameSubTree(TNode n1 , TNode n2){
        List<TNode> l = new ArrayList<TNode>();
        List<TNode> lSub = new ArrayList<TNode>();
        Iterator it = n1.getChildren().iterator();
        while(it.hasNext()){
            TNode temp = (TNode)it.next();
            if(temp.isLeaf())
                l.add(temp);
            else
                lSub.add(temp);
        }
        List<TNode> r = new ArrayList<TNode>();
        List<TNode> rSub = new ArrayList<TNode>();
        it = n2.getChildren().iterator();
        while(it.hasNext()){
            TNode temp = (TNode)it.next();
            if(temp.isLeaf())
                r.add(temp);
            else
                rSub.add(temp);
        }
        if(l.size()!=r.size())
            return false;
        else {
            if(l.size()==1){
                if(ifSameLeaf(l.get(0), r.get(0)))
                    return ifSameSubTree(lSub.get(0),rSub.get(0));
                else
                    return false;
            }
            else if(l.size()==2){
                if(ifSameLeaf(l.get(0), r.get(0)))
                    return ifSameLeaf(l.get(1), r.get(1));

                else if(ifSameLeaf(l.get(0), r.get(1)))
                    return ifSameLeaf(l.get(1), r.get(0));

                else
                    return false;
            }
            else{
                if(ifSameSubTree(lSub.get(0), rSub.get(0)))
                    return ifSameSubTree(lSub.get(1), rSub.get(1));

                else if(ifSameSubTree(lSub.get(0), rSub.get(1)))
                    return ifSameSubTree(lSub.get(1), rSub.get(0));

                else
                    return false;

            }
        }

    }

    public Tree getSubTree(Tree t, TNode n){
        TNode p = n.getParent();
        Tree temp = Trees.readTree(t.toNewick());
        temp.rerootTreeAtNode(n);
        return temp;
    }


    public List<String> getMSTopos(List<Tree> msGTS, List<String> checkedTopo){
        List<String> topos = new ArrayList<String>();
        for(int i = 0; i<msGTS.size(); i++){
            String temp = getTopo(msGTS.get(i).toString());
            if(checkedTopo.contains(temp))
                continue;
            else
                topos.add(getTopo(msGTS.get(i).toString()));
        }
        return topos;
    }

    public String getTopo(String tree){
        Tree temp = Trees.readTree(tree);
        Iterator node = temp.getNodes().iterator();
        while(node.hasNext()){
            TNode n = (TNode) node.next();
            if(!n.isRoot())
                //TODO 7 26
                n.setParentDistance(Double.NEGATIVE_INFINITY);
        }
        return temp.toString();
    }

    // ./raxmlHPC-PTHREADS -f e -t RAxML_bestTree.T0 -m GTRGAMMA -s 0/dna.phy -n TEST
   // rm *.TEST
//    public static List<Tree> addGTOG(List<Tree> gts, List<Double> ogHeight) {
//        List<Tree> newGTS = new ArrayList<Tree>();
//        for (int i = 0; i < gts.size(); i++) {
//            STITree thisTree = (STITree) gts.get(i);
//            STINode root = thisTree.getRoot();
////            STINode newroot = root.createChild("I" + thisTree.getNodeCount());
////            root.removeChild((TMutableNode) newroot,false);
////            newroot.adoptChild((TMutableNode)root);
////            STINode outgroup = newroot.createChild("O");
////            outgroup.setParentDistance(ogHeight.get(i));
////            root.setParentDistance(Double.NEGATIVE_INFINITY);
////            thisTree.rerootTreeAtNode(newroot);
////            newGTS.add((Tree)thisTree);
//
//            if (thisTree.getNode("O") != null) {
//                newGTS.add((Tree) thisTree);
//                continue;
//            }
//            STINode outgroup = root.createChild("O");
//            outgroup.setParentDistance(ogHeight.get(i));
//            newGTS.add((Tree) thisTree);
//        }
//        return newGTS;
//    }

//    public static List<Tree> removeGTOG(List<Tree> gts) {
//        List<Tree> newGTS = new ArrayList<Tree>();
//        for (int i = 0; i < gts.size(); i++) {
//            MutableTree thisTree = (MutableTree) gts.get(i);
//            TMutableNode root = thisTree.getRoot();
////            TNode outGroup = thisTree.getNode("O");
////            Iterator it = root.getChildren().iterator();
////            TNode newroot = (TNode)it.next();
////            if(newroot.getName().equals("O"))
////                newroot = (TNode)it.next();
////
////            thisTree.rerootTreeAtNode(newroot);
//            root.removeChild(thisTree.getNode("O"), false);
//            newGTS.add(thisTree);
//        }
//        return newGTS;
//    }

    public List<Tree> getOGGTS(List<Tree> gts, double ogHeight){
        List<Tree> ogGTS = new ArrayList<Tree>();
        for(int i = 0;i<gts.size();i++){
            String temp = gts.get(i).toString();
            int length = temp.length();
            Tree ogTree = Trees.readTree("(" + temp.substring(0,length-1) + ",O:" + String.valueOf(ogHeight) + ");");
            ogGTS.add(ogTree);
        }
        return ogGTS;
    }

    //DO NOT change original trees
    public List<Tree> getNoOGGTS(List<Tree> ogGTS){
        List<Tree> gts = new ArrayList<Tree>();
        for(int i = 0;i<ogGTS.size();i++){
            String temp = ogGTS.get(i).toString();
            Tree gt = Trees.readTree(rerootAndRemove(temp, "O"));
            gts.add(gt);
        }
        return gts;
    }

    //Input: rooted binary tree with oG in the first level
    //Output: rooted tree without OG
    public static String rerootAndRemove(String tree, String outgroup) {
        STITree t = (STITree) Trees.readTree(tree);
        if (t == null || t.getNode(outgroup) == null) {
            System.err.println(tree + "  " + outgroup);
            return null;
        }
        t.rerootTreeAtEdge(outgroup);
        t.removeNode(outgroup);
        Trees.removeBinaryNodes(t);
        return t.toNewick();
    }

    //reroot RAxML tree whose root has three children. DO NOT remove OG
    //Input: RAxML tree
    //Output: rooted tree and outgroup is the first-level child
    public static Tree rerootRAxML(Tree t, String outgroup) {
        Trees.autoLabelNodes((MutableTree) t);
        TNode og = t.getNode(outgroup);
        double oHeight = og.getParentDistance();
        String tempP = og.getParent().getName();
        TNode oldRoot = t.getRoot();
        int count = t.getNodeCount();
        if(oldRoot.getChildCount()!=3)
            System.err.println("Special RAxML tree!(2-node root)");
        BniNetwork tempNet = (BniNetwork)Networks.readNetwork(t.toString());
        BniNetNode newRoot = new BniNetNode("I" + String.valueOf(count),0.0);
        BniNetNode oldP = (BniNetNode)tempNet.findNode(tempP);

        oldP.adoptChild(newRoot,oHeight/100);
        newRoot.adoptChild((BniNetNode)tempNet.findNode(outgroup),oHeight);
        oldP.removeChild((BniNetNode)tempNet.findNode(outgroup));

        BniNetNode oldNetRoot = (BniNetNode)tempNet.getRoot();
        BniNetNode tempChild = newRoot;
        BniNetNode tempParent = oldP;
        while(!tempChild.equals(oldNetRoot)){
            Double tempD = tempChild.getParentDistance(tempParent);
            tempParent.removeChild(tempChild);
            tempChild.adoptChild(tempParent,tempD);
            tempChild = tempParent;
            tempParent = (BniNetNode) tempChild.getParents().iterator().next();
        }
        tempNet.resetRoot(newRoot);
        //Tree tempT = tempNet.
        t = Trees.readTree(tempNet.toString());
        return Trees.readTree(tempNet.toString());

    }

    //DO REMOVE OG
    public static Network reRootAST(Tree t, String outgroup) {

        Trees.autoLabelNodes((MutableTree)t);
        TMutableNode oldRoot = (TMutableNode)t.getRoot();
        if(oldRoot.getChildCount()!=2){
            System.err.print("Wrong child number for AST Tree!");
        }
        Iterator it = oldRoot.getChildren().iterator();
        TMutableNode newRoot = oldRoot;
        TMutableNode leafNode = oldRoot;
        while(it.hasNext()){
            TMutableNode temp = (TMutableNode)it.next();
            if(temp.isLeaf())
                leafNode = temp;
            else
                newRoot = temp;
        }
        newRoot.adoptChild(leafNode);
        t.rerootTreeAtNode(newRoot);
        newRoot.removeChild(oldRoot,false);
        Network st = Networks.readNetwork(rerootRAxML(t,outgroup).toNewick());
        NetNode og = st.findNode(outgroup);
        Iterator sibIt = st.getRoot().getChildren().iterator();
        NetNode sib = (NetNode) sibIt.next();
        if(sib.equals(og))
            sib = (NetNode) sibIt.next();
        st.resetRoot(sib);
        return st;
    }



    //UNDO: true -> *halfTheta    ;   false -> /halfTheta
    public static Tree scaleGT(Tree t, double ratio, boolean multiply) {

            Iterator itNode = t.getNodes().iterator();
            if (multiply == false) {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) t.getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() / ratio); // height = height * ratio
                    }
                }
            } else {
                while (itNode.hasNext()) {
                    STINode n = (STINode) itNode.next();
                    if (n == (STINode) t.getRoot())
                        continue;
                    if (n.getParentDistance() != TMutableNode.NO_DISTANCE) {
                        n.setParentDistance(n.getParentDistance() * ratio); // height = height / ratio
                    }
                }
            }
        return t;
    }

    public List<Tree> getTreeByUTree(List<UltrametricTree> uList) {
        Iterator uts = uList.iterator();
        List<Tree> gts = new ArrayList<Tree>();
        while (uts.hasNext()) {
            STITree gTree = new STITree(((UltrametricTree) uts.next()).getTree());
            gts.add(gTree);
        }
        return gts;
    }

    public List<UltrametricTree> getUTreeByTree(List<Tree> gt) {
        //List<UltrametricTree> ugts = new ArrayList<Tree>();
        Iterator gts = gt.iterator();
        List<UltrametricTree> ultraGTS = new ArrayList<UltrametricTree>();
        while (gts.hasNext()) {
            UltrametricTree uTree = new UltrametricTree((Tree) gts.next());
            ultraGTS.add(uTree);
        }
        return ultraGTS;
    }

    //Input: gt, shortest external length(ST height)
    public Tree getUTree(Tree gt, double base) {
        TNode lowLeaf = gt.getNodes().iterator().next();
        double distance = 0.0;
        Iterator nodeIt = gt.getNodes().iterator();
        while(nodeIt.hasNext()){
            TNode temp = (TNode)nodeIt.next();
            if(temp.isLeaf()){
                double tempDistance = getLeafHeight(temp, gt);
                if(tempDistance>distance){
                    distance = tempDistance;
                    lowLeaf = temp;
                }
            }
            else
                continue;
        }
        double height = distance + base;

        Iterator nodeIt1 = gt.getNodes().iterator();
        while(nodeIt1.hasNext()){
            TNode temp = (TNode)nodeIt1.next();
            if(temp.isLeaf()){
                double addition = height - getLeafHeight(temp, gt);
                temp.setParentDistance(temp.getParentDistance()+addition);

            }
        }
        //ltrametricTree ultraGT = new UltrametricTree(gt);
        return gt;
    }

    public Tree cleanName(Tree t) {
        Iterator node = t.getNodes().iterator();
        while(node.hasNext()){
            STINode stNode = (STINode)node.next();
            if(stNode.isLeaf())
                continue;
            stNode.setName("");
        }
        return t;
    }

    public Network getUTree(Network tempST, double[] height) {
        Tree st = Trees.readTree(tempST.toString());
        TNode lowLeaf = st.getNodes().iterator().next();
        double distance = 0.0;
        Iterator nodeIt = st.getNodes().iterator();
        while(nodeIt.hasNext()){
            TNode temp = (TNode)nodeIt.next();
            if(temp.isLeaf()){
                double tempDistance = getLeafHeight(temp, st);
                if(tempDistance>distance){
                    distance = tempDistance;
                    lowLeaf = temp;
                }
            }
            else
                continue;
        }
        height[2] = distance;

        Iterator nodeIt1 = st.getNodes().iterator();
        while(nodeIt1.hasNext()){
            TNode temp = (TNode)nodeIt1.next();
            if(temp.isLeaf()){
                double addition = height[2] - getLeafHeight(temp, st);
                temp.setParentDistance(temp.getParentDistance()+addition);

            }
        }
        return Networks.readNetwork(st.toString());
    }


    public double getDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, false);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double getRootDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, true);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    public double sumDoubleList(List<Double> list) {
        double sum = 0.0;
        Iterator it = list.iterator();
        while (it.hasNext()) {
            sum += (double) it.next();
        }
        return sum;
    }

    public double getLeafHeight(NetNode leaf, Network currentST) {
        double currentHeight = 0.0;
        //NetNode leaf = (NetNode) currentST.getLeaves().iterator().next();
        while (leaf != currentST.getRoot()) {
            currentHeight += leaf.getParentDistance((NetNode) leaf.getParents().iterator().next());
            leaf = (NetNode) leaf.getParents().iterator().next();
        }
        return currentHeight;
    }

    //for Tree
    public double getLeafHeight(TNode leaf, Tree currentST) {
        double currentHeight = 0.0;
        //NetNode leaf = (NetNode) currentST.getLeaves().iterator().next();
        while (leaf != currentST.getRoot()) {
            currentHeight += leaf.getParentDistance();
            leaf = (TNode) leaf.getParent();
        }
        return currentHeight;
    }

    public double getNodeHeight(Network currentST, double treeHeight, NetNode<Double> n) {
        //double currentHeight = getNetHeight(currentST);
        double height = 0.0;
        NetNode<Double> leaf = (NetNode<Double>) n.getParents().iterator().next();
        while (leaf != currentST.getRoot()) {
            height += leaf.getParentDistance((NetNode<Double>) leaf.getParents().iterator().next());
            leaf = (NetNode<Double>) leaf.getParents().iterator().next();
        }
        n.setParentDistance((NetNode<Double>) n.getParents().iterator().next(),treeHeight-height);
        return treeHeight-height;
    }

    //
    //
    ////////////////////////////////////
    //For UPRA:

    //Input: MS sample trees
    //Output: distinguished MS topos
    public List<Tree> getDistinguishedTopo(List<Tree> msTrees){
        List<Tree> dTopos = new ArrayList<Tree>();
        List<Integer> sameTopoNum = new ArrayList<Integer>();

        for(int i = 0; i < msTrees.size(); i++){
            if(sameTopoNum.contains(i))
                continue;
            Tree t1 = msTrees.get(i);
            dTopos.add(t1);

            for(int j = i+1 ; j < msTrees.size(); j++){
                if(sameTopoNum.contains(j))
                    continue;
                Tree t2 = msTrees.get(j);
                if(Trees.haveSameRootedTopology(t1,t2)){
                    sameTopoNum.add(j);
                }
            }

        }

        return dTopos;
    }



}
