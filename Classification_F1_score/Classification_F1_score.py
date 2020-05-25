import laspy
import numpy as np
from sklearn.metrics import f1_score, confusion_matrix, classification_report
from tabulate import tabulate
from multiprocessing import Pool


def get_labels(file): 
    pc = laspy.file.File(file, mode='r')
    ids = pc.gps_time
    return pc.classification, list(ids)

def get_labels_CSF(file): 
    pc = laspy.file.File(file, mode='r')
    ids = []
    for i in range(pc.__len__()):
      X = pc.X[i] * pc.header.scale[0] + pc.header.offset[0]
      Y = pc.Y[i] * pc.header.scale[1] + pc.header.offset[1]
      Z = pc.Z[i] * pc.header.scale[2] + pc.header.offset[2]
      ids.append([X, Y, Z])
    return pc.classification, list(ids)


def f1score(true_labels, predicted_labels, file):
    obj = open(file, 'w')
    
    for i in range(len(true_labels)):
        if true_labels[i] != 2:
            true_labels[i] = 0
    for i in range(len(predicted_labels)):
        if predicted_labels[i] != 2:
            predicted_labels[i] = 0

    con_matrix = confusion_matrix(true_labels, predicted_labels,[2,0])
    class_report = classification_report(true_labels, predicted_labels)

    table = [ ["", "Actual ground points", "Actual non-ground points"], 
            ["Predicted ground points", con_matrix[0][0], con_matrix[0][1]], 
            ["Predicted non-ground points", con_matrix[1][0], con_matrix[1][1]] ]
    
    s1 = f1_score(true_labels, predicted_labels, average=None)[1]
    s2 = f1_score(true_labels, predicted_labels, average='macro')
    s3 = f1_score(true_labels, predicted_labels, average='micro')
    s4 = f1_score(true_labels, predicted_labels, average='weighted')
    
    obj.write("Confusion matrix\n")
    obj.write(tabulate(table))
    obj.write("\n\nF1-score = " + str(s1) + "\n")
    obj.write("macro F1-score = " + str(s2) + "\n")
    obj.write("micro F1-score = " + str(s3) + "\n")
    obj.write("weighted F1-score = " + str(s4))
    obj.write("\n\nOverall results table\n")
    obj.write("-------------------------------------\n\n")
    obj.write(class_report)

    obj.close()
    return

def find_points( parameters ):
    time_true = parameters[0]
    time_pre = parameters[1]
    labels_true = parameters[2]
    labels_pre = parameters[3]
    new_true_labels=[]
    new_predicted_labels=[]
    for i in range(len(time_true)):
        if time_true[i] in time_pre:
            new_true_labels.append( labels_true[i] )
            new_predicted_labels.append( labels_pre[time_pre.index(time_true[i])] )
    
    return [new_true_labels, new_predicted_labels]

def find_and_fix_points( parameters ):
    time_true = parameters[0]
    time_pre = parameters[1]
    labels_true = parameters[2]
    new_true_labels=[]
    new_predicted_labels=[]

    for i in range(len(time_true)):
        if time_true[i] in time_pre:
            new_true_labels.append( labels_true[i] )
            new_predicted_labels.append( 2 )
        else:
            new_true_labels.append( labels_true[i] )
            new_predicted_labels.append( 0 )
    
    return [new_true_labels, new_predicted_labels]


if __name__ == '__main__':
    true_labels_ams, a = get_labels("true_labels/Amsterdam_manually_classified.las")
    true_labels_del, a2 = get_labels("true_labels/Delft_manually_classified.las")
    true_labels_bie, a3 = get_labels("true_labels/Biesbosch_manually_classified.las")
    gf_methods = ["_AHN", "_lasground"]
    
    for method in gf_methods:
        predicted_labels_ams, b = get_labels(method[1:] + "/Amsterdam" + method + ".las")
        predicted_labels_del, b2 = get_labels(method[1:] + "/Delft" + method + ".las")
        predicted_labels_bie, b3 = get_labels(method[1:] + "/Biesbosch" + method + ".las")

        f1score(true_labels_ams, predicted_labels_ams, method[1:] + "/Amsterdam" + method + "_report.txt")
        f1score(true_labels_del, predicted_labels_del, method[1:] + "/Delft" + method + "_report.txt")
        f1score(true_labels_bie, predicted_labels_bie, method[1:] + "/Biesbosch" + method + "_report.txt")

    problematic_gf_methods = ["_PDAL_smrf", "_PDAL_pmf"]
    for method in problematic_gf_methods:
        predicted_labels_ams, c = get_labels(method[1:] + "/Amsterdam" + method + ".las")
        predicted_labels_del, c2 = get_labels(method[1:] + "/Delft" + method + ".las")
        predicted_labels_bie, c3 = get_labels(method[1:] + "/Biesbosch" + method + ".las")
        pre_map = []
        pre_map.append([a, c, true_labels_ams, predicted_labels_ams])
        pre_map.append([a2, c2, true_labels_del, predicted_labels_del])
        pre_map.append([a3, c3, true_labels_bie, predicted_labels_bie])
        p = Pool(processes = 3)
        labels = p.map(find_points, pre_map)
        p.close(); p.join() 
        f1score(labels[0][0], labels[0][1], method[1:] + "/Amsterdam" + method + "_report.txt")
        f1score(labels[1][0], labels[1][1], method[1:] + "/Delft" + method + "_report.txt")
        f1score(labels[2][0], labels[2][1], method[1:] + "/Biesbosch" + method + "_report.txt")
    
    true_labels_ams, a = get_labels_CSF("true_labels/Amsterdam_manually_classified.las")
    true_labels_del, a2 = get_labels_CSF("true_labels/Delft_manually_classified.las")
    true_labels_bie, a3 = get_labels_CSF("true_labels/Biesbosch_manually_classified.las")
    i_do_whatever_I_want_method = ["_CSF"] 
    # gives as a result an las file with only the ground points and does not change their class to 2, the points have multiple classes
    # also deletes the gps_time which can be used as id for every point so the coordinates had to be used
    for method in i_do_whatever_I_want_method:
        predicted_labels_ams, d = get_labels_CSF(method[1:] + "/Amsterdam" + method + ".las")
        predicted_labels_del, d2 = get_labels_CSF(method[1:] + "/Delft" + method + ".las")
        predicted_labels_bie, d3 = get_labels_CSF(method[1:] + "/Biesbosch" + method + ".las")
        pre_map = []
        pre_map.append([a, d, true_labels_ams])
        pre_map.append([a2, d2, true_labels_del])
        pre_map.append([a3, d3, true_labels_bie])
        p = Pool(processes = 3)
        labels = p.map(find_and_fix_points, pre_map)
        p.close(); p.join() 
        f1score(labels[0][0], labels[0][1], method[1:] + "/Amsterdam" + method + "_report.txt")
        f1score(labels[1][0], labels[1][1], method[1:] + "/Delft" + method + "_report.txt")
        f1score(labels[2][0], labels[2][1], method[1:] + "/Biesbosch" + method + "_report.txt")