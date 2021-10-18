## Import Lines
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn import model_selection as ms

##Create empty variables
target = []
not_target = []
unclassified = []


## reads in data from the C scripts. 
## seqname = contig name, seqlen = contig length, gc = GC content, cov = coverage, hasBlast = integer for Blast hit presence, isTarget = integer for Blast hit being target or not-target, kdist_list= kmer list
## This function sorts out the data into target, not_target, and unclassified. 
def import_seqdata(seqname, seqlen, gc, cov, hasBlast, isTarget, kdist_list):
    seqdata = (seqname, seqlen, kdist_list, isTarget, gc, cov)
#    print("Seqdata imported %d " % len(seqdata))
    if hasBlast == 1:
        if isTarget == 1:
            target.append(seqdata)
        else:
            not_target.append(seqdata)
    else:
        unclassified.append(seqdata)
    return 0



def DT_classify(classifier):
    regionIDs = []
    X = []
    for item in unclassified:
        regionIDs.append(item[0])       # region name
        X.append(item[4:])
    Y = []
    for i in classifier.predict(X):
        Y.append(i)
    return list(zip(regionIDs, Y))


## Decision Tree Function. Uses target and not_target to train. Then Classifies unclassified contigs
def DT_model(unclassified):
    X = []
    Y = []
    features = ["GC", "Coverage"]
    for item in target + not_target:
        X.append(item[4:])
        Y.append(item[3])       # isTarget
    X_train, X_test, Y_train, Y_test = ms.train_test_split(X, Y, test_size=0.33, random_state=0)
    classifier = tree.DecisionTreeClassifier()
    classifier = classifier.fit(X_train, Y_train)
    print("Classifier built, score is %s out of 1.00" % classifier.score(X_test, Y_test))
    #with open("model.dot", 'w') as dotfile:
    #    tree.export_graphviz(classifier, out_file=dotfile, feature_names=features,
    #                         class_names=Y, filled=True, rounded=True, special_characters=True)
    return DT_classify(classifier)


## Print statements for out, calls decision Tree function, writes final output
def run_analysis(kmer_list):
    print("Seqdata imported, %d target sequences, %d contaminant sequences, and %d unclassified sequences\n" % (len(target), len(not_target), len(unclassified)))

## Call Decision Tree
    if len(unclassified) != 0:
        result = DT_model(unclassified)
    else:
        result = []
        print("No unclassified sequences to analyze")

## write tokeep.txt, toremove.txt, and fulldata.txt
## tokeep and toremove are lists of sequences and their label. tab delinated
## full data is all of the sequences, their labels, GC Content, and Coverage. Tab deliniated
    with open("tokeep.txt", "w+") as k, open("toremove.txt", "w+") as r, open("fulldata.txt", "w+") as f:
        k.write("Sequence_Name\tLabel\n") #header
        r.write("Sequence_Name\tLabel\n") #header
        f.write("Sequence_Name\tLabel\tGC_Content\tCoverage\n") #header 
        for item in target:
           k.write("%s\ttarget\n" % item[0])
           f.write("%s\ttarget\t%s\t%s\n" % (item[0], str(item[4]), str(item[5])))
        for item in not_target:
           r.write("%s\tcontaminant\n" % item[0])
           f.write("%s\tcontaminant\t%s\t%s\n" % (item[0], str(item[4]), str(item[5])))
        for item in result:
            if item[1] == 1:
                k.write("%s\tclassified_target\n" % item[0])
                print("%s was classified as Target\n" % item[0])
                for line in unclassified:
                    if item[0] == line[0]:
                        f.write("%s\tclassified_target\t%s\t%s\n" % (line[0], str(line[4]), str(line[5])))
            else:
                r.write("%s\tclassified_contaminant\n" % item[0])
                print("%s was classified as Contaminant\n" % item[0])
                for line in unclassified:
                    if item[0] == line[0]:
                        f.write("%s\tclassified_contaminant\t%s\t%s\n" % (line[0], str(line[4]), str(line[5])))
    return 1
