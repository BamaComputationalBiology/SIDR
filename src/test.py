from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn import model_selection as ms



target = []
not_target = []
unclassified = []


def import_seqdata(seqname, seqlen, gc, cov, hasBlast, isTarget, kdist_list):
    seqdata = (seqname, seqlen, kdist_list, isTarget, gc, cov)
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



def run_analysis(kmer_list):
    print("Seqdata imported, %d target sequences, %d contaminant sequences, and %d unclassified sequences\n" % (len(target), len(not_target), len(unclassified)))
    if len(unclassified) != 0:
        result = DT_model(unclassified)
    else:
        result = []
    with open("tokeep.txt", "w+") as k, open("toremove.txt", "w+") as r:
        for item in result:
            if item[1] == 1:
                k.write("%s x\n" % item[0])
            else:
                r.write("%s x\n" % item[0])

    with open("target.txt", "w+") as t, open("nottarget.txt", "w+") as n:
	for item in target:
		t.write("%s x\n" % item)
	for item in nottarget:
		n.write("%s x\n" % item)


    return 1
