from PyQt5 import QtWidgets, QtGui
from PyQt5.QtWidgets import *
import sys
import os
from ip_processing import start_pool


class MyWindow(QMainWindow):
    def __init__(self):
        super(MyWindow, self).__init__()
        self.initUI()

    def select_directory(self):
        selected_dir = QFileDialog.getExistingDirectory(self, caption='Choose Directory', directory=os.getcwd())
        self.folder.setText(selected_dir+"/")
        return

    def button_clicked(self):
        lst = ["no post-processing", "patching","polygon flattening", 'patching & polygon flattening']
        for i in lst:
            if self.postprocess.currentText() == i:
                print(int(lst.index(i)))
                postprocess_parameter = int(lst.index(i))
        idw0_p = 5
        idw1_p = 2
        idw2_p = 0
        idw3_p = 2
        idw5_p = 0.2
        idw6_p = 3
        size_p = 1

        if self.idw0.text() != "":
            idw0_p = int(self.idw0.text())

        if self.idw1.text() != "":
            idw1_p = int(self.idw1.text())

        if self.idw2.text() != "":
            idw2_p = int(self.idw2.text())

        if self.idw3.text() != "":
            idw3_p = int(self.idw3.text())

        if self.idw5.text() != "":
            idw5_p = int(self.idw5.text())

        if self.idw6.text() != "":
            idw6_p = int(self.idw6.text())

        if self.idw5.text() != "":
            size_p = float(self.idw5.text())

        start_pool(self.folder.text(), preprocess = self.preprocess.currentText(), postprocess = postprocess_parameter,
                size = size_p, method = self.method.currentText(), fmt = self.fmt.currentText(),
                idw0_polyfpath = idw0_p, idw1 = idw1_p, idw2 = idw2_p, idw3 = idw3_p,
                idw4 = self.idw4.currentText(), idw5 = idw5_p, idw6 = idw6_p)

        msg = QMessageBox()
        msg.setWindowTitle("Information")
        msg.setText("Process done!")
        msg.setIcon(QMessageBox.Information)
        msg.setStandardButtons(QMessageBox.Ok)
        text = "Parameters used:\n"+ self.method.currentText() + "\n" + self.preprocess.currentText() + "\n" + self.postprocess.currentText() + "\n" + self.fmt.currentText() + "\n "+ str(size_p) + "\n" + str(idw0_p) + "\n"+ str(idw1_p) + "\n"+ str(idw2_p) + "\n" + str(idw3_p) + "\n"+ self.idw4.currentText() + "\n" + str(idw5_p) + "\n" + str(idw6_p)
        msg.setDetailedText(text)
        msg.exec_()
        return

    def initUI(self):
        
        self.setGeometry(100, 100, 800, 800)
        self.setWindowTitle("Interpolation GUI")

        self.label = QtWidgets.QLabel(self)
        self.label.setText("Interpolation parameters")
        self.label.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold))
        self.label.move(250, 20)
        self.label.adjustSize()
        
        self.folder_label = QtWidgets.QLabel(self)
        self.folder_label.setText("Select the target folder:")
        self.folder_label.move(100, 100)
        self.folder_label.adjustSize()
        self.folder = QLineEdit(self)
        self.folder.move(500, 100)
        self.folder.resize(140, 30)
        
        self.b0 = QtWidgets.QPushButton(self)
        self.b0.setText("...")
        self.b0.clicked.connect(self.select_directory)
        self.b0.move(640, 100)
        self.b0.resize(35, 30)

        self.method_label = QtWidgets.QLabel(self)
        self.method_label.setText("Choose an interpolation method:")
        self.method_label.move(100, 150)
        self.method_label.adjustSize()
        self.method = QComboBox(self)
        self.method.addItems(["startin-Laplace","PDAL-IDW", "startin-TINlinear", 'CGAL-NN', 'CGAL-CDT', 'IDWquad'])
        self.method.move(500, 150)
        self.method.resize(175, 30)

        self.preprocess_label = QtWidgets.QLabel(self)
        self.preprocess_label.setText("Run a ground filtering algorithm?")
        self.preprocess_label.move(100, 200)
        self.preprocess_label.adjustSize()
        self.preprocess = QComboBox(self)
        self.preprocess.addItems(["False", "True"])
        self.preprocess.move(500, 200)
        self.preprocess.resize(175, 30)

        self.postprocess_label = QtWidgets.QLabel(self)
        self.postprocess_label.setText("What post-processing method implement?")
        self.postprocess_label.move(100, 250)
        self.postprocess_label.adjustSize()
        self.postprocess = QComboBox(self)
        self.postprocess.addItems(["no post-processing", "patching","polygon flattening", 'patching & polygon flattening'])
        self.postprocess.move(500, 250)
        self.postprocess.resize(175, 30)

        self.fmt_label = QtWidgets.QLabel(self)
        self.fmt_label.setText("Export format:")
        self.fmt_label.move(100, 300)
        self.fmt_label.adjustSize()
        self.fmt = QComboBox(self)
        self.fmt.addItems(["GeoTIFF", "ASC"])
        self.fmt.move(500, 300)
        self.fmt.resize(175, 30)

        self.size_label = QtWidgets.QLabel(self)
        self.size_label.setText("Raster cell size (meters):")
        self.size_label.move(100, 350)
        self.size_label.adjustSize()
        self.size = QLineEdit(self)
        self.size.move(500, 350)
        self.size.resize(175,30)

        self.idw0_label = QtWidgets.QLabel(self)
        self.idw0_label.setText("IDW radius (meters):")
        self.idw0_label.move(100, 400)
        self.idw0_label.adjustSize()
        self.idw0 = QLineEdit(self)
        self.idw0.move(500, 400)
        self.idw0.resize(175,30)

        self.idw1_label = QtWidgets.QLabel(self)
        self.idw1_label.setText("IDW power:")
        self.idw1_label.move(100, 450)
        self.idw1_label.adjustSize()
        self.idw1 = QLineEdit(self)
        self.idw1.move(500, 450)
        self.idw1.resize(175,30)

        self.idw2_label = QtWidgets.QLabel(self)
        self.idw2_label.setText("IDW fallback kernel width (for PDAL-IDW):")
        self.idw2_label.move(100, 500)
        self.idw2_label.adjustSize()
        self.idw2 = QLineEdit(self)
        self.idw2.move(500, 500)
        self.idw2.resize(175,30)

        self.idw3_label = QtWidgets.QLabel(self)
        self.idw3_label.setText("Radius/number of neighbours INCREMENT value (for IDWquad):")
        self.idw3_label.move(100, 550)
        self.idw3_label.adjustSize()
        self.idw3 = QLineEdit(self)
        self.idw3.move(500, 550)
        self.idw3.resize(175,30)

        self.idw4_label = QtWidgets.QLabel(self)
        self.idw4_label.setText("IDWquad method:")
        self.idw4_label.move(100, 600)
        self.idw4_label.adjustSize()
        self.idw4 = QComboBox(self)
        self.idw4.addItems(["radial", "k-nearest"])
        self.idw4.move(500, 600)
        self.idw4.resize(175,30)

        self.idw5_label = QtWidgets.QLabel(self)
        self.idw5_label.setText("IDWquad tolerance value (epsilon):")
        self.idw5_label.move(100, 650)
        self.idw5_label.adjustSize()
        self.idw5 = QLineEdit(self)
        self.idw5.move(500, 650)
        self.idw5.resize(175,30)

        self.idw6_label = QtWidgets.QLabel(self)
        self.idw6_label.setText("IDWquad maximum number of iteration before declaring no-data:")
        self.idw6_label.move(100, 700)
        self.idw6_label.adjustSize()
        self.idw6 = QLineEdit(self)
        self.idw6.move(500, 700)
        self.idw6.resize(175,30)

        self.b1 = QtWidgets.QPushButton(self)
        self.b1.setText("Start interpolation")
        self.b1.clicked.connect(self.button_clicked)
        self.b1.move(250, 750)
        self.b1.resize(150,30)
     
        self.show()


def window():
    app = QApplication(sys.argv)
    win = MyWindow()
    win.show()
    sys.exit(app.exec_())



if __name__ == '__main__':
    window()