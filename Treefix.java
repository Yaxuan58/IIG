package edu.rice.cs.bioinfo.programs.phylonet.algos.iterHeuristic;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by doriswang on 1/23/18.
 */
public class Treefix {

    private String _TreefixDir;
    private String _InputDir;
    private String _OutputDir;
    private String _RAxMLDir;

    private String _RTreeDir; // rooted raxml output
    private int TREEINDEX;
    private static List<Tree> trueGTS;
    private static Tree trueST;
    private static List<Alignment> trueSeq;
    private int lociNum;
    private InferOperator operator;
    private int seqLen;
    private int taxaNum;
    private String theta;
    private static String gtPath;


    Treefix() throws IOException, ParseException, InterruptedException {


        TREEINDEX = 0;
        trueGTS = new ArrayList<Tree>();
        trueSeq = new ArrayList<Alignment>();

        //change variables
        lociNum = 200;
        seqLen = 1000;
        taxaNum = 16;
        theta = "001";


        gtPath = "/Users/doriswang/PhyloNet/tools/RAxML/" + seqLen + "t" + theta + "/";

        _TreefixDir = "/Users/doriswang/Treefix/venv/treefix-1.1.10/s16/";
        _InputDir = "/Users/doriswang/PhyloNet/";
        _OutputDir = "/Users/doriswang/PhyloNet/output/RAxML/" + theta + "/"; //for RAxML out
        _RTreeDir = _OutputDir + "rooted/";

        operator = new InferOperator(_InputDir, _OutputDir, 0, 7);
        operator.isExitsPath(_OutputDir);
        _RAxMLDir = operator._RAxMLdir;
        loadTrees(_InputDir + "input/" + theta + "/", lociNum);
        for (int i = 0; i < lociNum; i++) {
            Alignment aln = operator.loadLocus(i, seqLen, taxaNum, _InputDir + "input/" + theta + "/");
            trueSeq.add(aln);
        }

    }

    //(1) Write alignment + update .stree
    //(2) Run Treefix in shell(tf.sh)
    //(3) Compute distance

    public static void main(String[] args) throws IOException, ParseException, InterruptedException {

        Treefix tf = new Treefix();
        tf.writeTFAlgn();


        String dataPath = gtPath;

        List<Tree> rTrees = loadRTrees(dataPath, tf.lociNum);

        double[] distR = tf.getDistances(trueGTS, rTrees, tf.lociNum);

        tf.writeRTrees(rTrees, tf._TreefixDir);

        List<Tree> tTrees = loadTFTrees(tf._TreefixDir, tf.lociNum);

        double[] distT = tf.getDistances(trueGTS, tTrees, tf.lociNum);

        BufferedWriter w = new BufferedWriter(new FileWriter(dataPath + "distanceR.txt"));
        BufferedWriter w1 = new BufferedWriter(new FileWriter(dataPath + "distanceTF.txt"));
        for (int i = 0; i < distR.length; i++) {
            w.write(i + ":" + distR[i] + "\n");
            w1.write(i + ":" + distT[i] + "\n");
            w.flush();
            w1.flush();
        }

        w.close();
        w1.close();
        System.out.println("R_Distance:" + distR[tf.lociNum]);
        System.out.println("TF_Distance:" + distT[tf.lociNum]);


        double distRS = 0.0;
        double distTS = 0.0;

        BufferedWriter w2 = new BufferedWriter(new FileWriter(dataPath + "distanceRS.txt"));
        BufferedWriter w3 = new BufferedWriter(new FileWriter(dataPath + "distanceTFS.txt"));
        for (int i = 0; i < distR.length - 1; i++) {
            distRS += tf.getDistance(trueST, rTrees.get(i));
            distTS += tf.getDistance(trueST, tTrees.get(i));
            w2.write(i + ":" + tf.getDistance(trueST, rTrees.get(i)) + "\n");
            w3.write(i + ":" + tf.getDistance(trueST, tTrees.get(i)) + "\n");
            w2.flush();
            w3.flush();
        }
        distRS = distRS / tf.lociNum;
        distTS = distTS / tf.lociNum;
        w2.write(tf.lociNum + ":" + distRS + "\n");
        w3.write(tf.lociNum + ":" + distTS + "\n");
        w2.flush();
        w3.flush();
        w2.close();
        w3.close();

        System.out.println("RS_Distance: " + distRS);
        System.out.println("TFS_Distance: " + distTS);


        List<Tree> t = loadTFTrees(tf._TreefixDir, tf.lociNum);
        w = new BufferedWriter(new FileWriter(tf._TreefixDir + tf.seqLen + "_" + tf.theta + "_TFTrees.txt"));
        w1 = new BufferedWriter(new FileWriter(tf._TreefixDir + tf.seqLen + "_" + tf.theta + "_RTrees.txt"));
        for (int i = 0; i < tf.lociNum; i++) {
            w.write(t.get(i).toString());
            w1.write(rTrees.get(i).toString());
            w.flush();
            w1.flush();
        }
        w.close();
        w1.close();
    }

    public void writeTFAlgn() throws IOException, ParseException, InterruptedException {
        for (int i = 0; i < lociNum; i++) {
            Alignment aln = trueSeq.get(i);
            String tempPath = _TreefixDir + i + "/";
            operator.isExitsPath(tempPath);
            BufferedWriter aFile = new BufferedWriter(new FileWriter(tempPath + "0.nt.align")); //0.nt.align
            List<String> names = aln.getTaxaNames();
            Map<String, String> thisAln = aln.getAlignment();
            for (int j = 0; j < taxaNum; j++) {
                String name = names.get(j);
                String seq = thisAln.get(name);
                aFile.write(">" + name + '\n');
                aFile.flush();
                aFile.write(seq + '\n');
                aFile.flush();
            }
            aFile.flush();
            aFile.close();
        }
    }

    public static List<Tree> loadTrees(String filePath, int lociNum) throws IOException {
        String streeFile = filePath + "Tree/trueST.txt";
        BufferedReader stReader = new BufferedReader(new FileReader(streeFile));
        trueST = Trees.readTree((String) stReader.readLine().trim());

        for (int ln = 0; ln < lociNum; ln++) {
            trueGTS.add(Trees.readTree(stReader.readLine().trim()));
        }
        return trueGTS;
    }

    //Input: filePath: /Users/doriswang/PhyloNet/tools/RAxML/TF200t3/

    public static List<Tree> loadRTrees(String filePath, int lociNum) throws IOException {
        String streeFile = filePath + "RAxML_bestTree.T";//RAxML_bestTree.T49
        ArrayList<Tree> rTrees = new ArrayList<Tree>();

        for (int ln = 0; ln < lociNum; ln++) {
            BufferedReader stReader = new BufferedReader(new FileReader(streeFile + ln));
            Tree temp = Trees.readTree(stReader.readLine().trim());
            rTrees.add(temp);
        }
        return rTrees;
    }

    //Users/doriswang/Treefix/venv/treefix-1.1.10/s16/2/0.nt.raxml.tree
    public static List<Tree> loadTFTrees(String filePath, int lociNum) throws IOException {

        ArrayList<Tree> tfTrees = new ArrayList<Tree>();

        for (int ln = 0; ln < lociNum; ln++) {

            String tree = readToString(filePath + ln + "/0.nt.raxml.treefix.tree");
            Tree temp = Trees.readTree(tree);
            tfTrees.add(temp);
        }
        return tfTrees;
    }

    public static String readToString(String filePath) {
        File file = new File(filePath);
        Long filelength = file.length();
        byte[] filecontent = new byte[filelength.intValue()];
        try {
            FileInputStream in = new FileInputStream(file);
            in.read(filecontent);
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new String(filecontent);
    }

    public double[] getDistances(List<Tree> l1, List<Tree> l2, int lociNum) throws IOException {
        Double dist = 0.0;
        double d[] = new double[lociNum + 1];

        for (int i = 0; i < l1.size(); i++) {
            d[i] = operator.getDistance(l1.get(i), l2.get(i));
            dist += operator.getDistance(l1.get(i), l2.get(i));
        }
        d[lociNum] = dist / lociNum;
        return d;
    }

    public double getDistance(Tree tree1, Tree tree2) {
        SymmetricDifference symmetricDifference = new SymmetricDifference();
        symmetricDifference.computeDifference(tree1, tree2, false);
        double diff = symmetricDifference.getWeightedAverage();
        return diff;
    }

    ///Users/doriswang/Treefix/venv/treefix-1.1.10/s16/2/0.nt.raxml.tree

    public void writeRTrees(List<Tree> rTrees, String tfPath) throws IOException {

        for (int i = 0; i < rTrees.size(); i++) {
            BufferedWriter w = new BufferedWriter(new FileWriter(tfPath + i + "/0.nt.raxml.tree"));
            Tree temp = rTrees.get(i);
            TNode node = (TNode) temp.getNodes().iterator().next();
            temp.rerootTreeAtNode(node);
            w.write(temp.toString());
            w.flush();
            w.close();
        }
    }
}



