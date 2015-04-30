"""
This module defines a class Ui_Form, consisiting of various attributes which corresponds to parameters of the cell. The
instance of this class is gives the User Interface which takes the input from the users and assigns it to the class
attributes.
In the interface, various SLOT-SIGNAL connects are used to make options exclusive. The function Process processes all the 
input. 


"""
from PyQt4 import QtCore, QtGui
import sys

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(QtGui.QWidget):

    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setupUi(self)
        self.time = 0.0
        self.anode_width = 0.0
        self.cathode_width = 0.0
        self.separator_width = 0.0
        self.conductivity_anode = 0.0
        self.porosity_anode = 0.0
        self.duffusion_coeff_anode = 0.0
        self.duffusion_coeff__solid_anode = 0.0
        self.reaction_constant_anode = 0.0
        self.conductivity_cathode = 0.0
        self.porosity_cathode = 0.0
        self.duffusion_coeff_cathode = 0.0
        self.duffusion_coeff__solid_cathode = 0.0
        self.reaction_constant_cathode = 0.0
        self.duffusion_coeff_separator = 0.0
        self.reaction_constant_separator = 0.0
        self.current_density = 0.0

    def setupUi(self, Form): 
	""" 
	This function sets up the User Interface and connects various SLOTs and SIGNALs
	"""
        Form.setObjectName(_fromUtf8("Form"))
        Form.setEnabled(True)
        Form.resize(1139, 713)
        Form.setAcceptDrops(False)
        self.label = QtGui.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(40, 50, 441, 31))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(40, 90, 83, 31))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.lineEdit = QtGui.QLineEdit(Form)
        self.lineEdit.setGeometry(QtCore.QRect(110, 90, 41, 31))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.label_3 = QtGui.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(270, 90, 83, 31))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.lineEdit_2 = QtGui.QLineEdit(Form)
        self.lineEdit_2.setGeometry(QtCore.QRect(360, 90, 41, 31))
        self.lineEdit_2.setObjectName(_fromUtf8("lineEdit_2"))
        self.lineEdit_3 = QtGui.QLineEdit(Form)
        self.lineEdit_3.setGeometry(QtCore.QRect(560, 90, 41, 31))
        self.lineEdit_3.setObjectName(_fromUtf8("lineEdit_3"))
        self.label_4 = QtGui.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(460, 90, 83, 31))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.checkBox = QtGui.QCheckBox(Form)
        self.checkBox.setGeometry(QtCore.QRect(710, 80, 191, 41))
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.line = QtGui.QFrame(Form)
        self.line.setGeometry(QtCore.QRect(0, 130, 981, 21))
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.checkBox_2 = QtGui.QCheckBox(Form)
        self.checkBox_2.setGeometry(QtCore.QRect(20, 140, 261, 71))
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.label_9 = QtGui.QLabel(Form)
        self.label_9.setEnabled(False)
        self.label_9.setGeometry(QtCore.QRect(30, 290, 111, 31))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.lineEdit_4 = QtGui.QLineEdit(Form)
        self.lineEdit_4.setEnabled(False)
        self.lineEdit_4.setGeometry(QtCore.QRect(150, 290, 41, 31))
        self.lineEdit_4.setObjectName(_fromUtf8("lineEdit_4"))
        self.label_10 = QtGui.QLabel(Form)
        self.label_10.setEnabled(False)
        self.label_10.setGeometry(QtCore.QRect(220, 290, 71, 31))
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.lineEdit_5 = QtGui.QLineEdit(Form)
        self.lineEdit_5.setEnabled(False)
        self.lineEdit_5.setGeometry(QtCore.QRect(310, 290, 41, 31))
        self.lineEdit_5.setObjectName(_fromUtf8("lineEdit_5"))
        self.label_11 = QtGui.QLabel(Form)
        self.label_11.setEnabled(False)
        self.label_11.setGeometry(QtCore.QRect(370, 290, 181, 31))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.lineEdit_6 = QtGui.QLineEdit(Form)
        self.lineEdit_6.setEnabled(False)
        self.lineEdit_6.setGeometry(QtCore.QRect(540, 290, 41, 31))
        self.lineEdit_6.setObjectName(_fromUtf8("lineEdit_6"))
        self.label_12 = QtGui.QLabel(Form)
        self.label_12.setEnabled(False)
        self.label_12.setGeometry(QtCore.QRect(610, 290, 181, 31))
        self.label_12.setAutoFillBackground(False)
        self.label_12.setWordWrap(True)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.lineEdit_7 = QtGui.QLineEdit(Form)
        self.lineEdit_7.setEnabled(False)
        self.lineEdit_7.setGeometry(QtCore.QRect(800, 290, 41, 31))
        self.lineEdit_7.setObjectName(_fromUtf8("lineEdit_7"))
        self.checkBox_4 = QtGui.QCheckBox(Form)
        self.checkBox_4.setGeometry(QtCore.QRect(330, 140, 261, 71))
        self.checkBox_4.setObjectName(_fromUtf8("checkBox_4"))
        self.line_3 = QtGui.QFrame(Form)
        self.line_3.setGeometry(QtCore.QRect(0, 520, 1111, 16))
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.label_5 = QtGui.QLabel(Form)
        self.label_5.setGeometry(QtCore.QRect(10, 550, 191, 31))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.lineEdit_8 = QtGui.QLineEdit(Form)
        self.lineEdit_8.setGeometry(QtCore.QRect(210, 550, 71, 31))
        self.lineEdit_8.setObjectName(_fromUtf8("lineEdit_8"))
        self.checkBox_5 = QtGui.QCheckBox(Form)
        self.checkBox_5.setGeometry(QtCore.QRect(450, 530, 361, 71))
        self.checkBox_5.setObjectName(_fromUtf8("checkBox_5"))
        self.lineEdit_9 = QtGui.QLineEdit(Form)
        self.lineEdit_9.setEnabled(True)
        self.lineEdit_9.setGeometry(QtCore.QRect(700, 550, 81, 31))
        self.lineEdit_9.setObjectName(_fromUtf8("lineEdit_9"))
        self.label_14 = QtGui.QLabel(Form)
        self.label_14.setGeometry(QtCore.QRect(290, 550, 191, 31))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.pushButton = QtGui.QPushButton(Form)
        self.pushButton.setGeometry(QtCore.QRect(470, 620, 151, 61))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.label_13 = QtGui.QLabel(Form)
        self.label_13.setEnabled(False)
        self.label_13.setGeometry(QtCore.QRect(30, 250, 111, 31))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.lineEdit_10 = QtGui.QLineEdit(Form)
        self.lineEdit_10.setEnabled(False)
        self.lineEdit_10.setGeometry(QtCore.QRect(1030, 290, 41, 31))
        self.lineEdit_10.setObjectName(_fromUtf8("lineEdit_10"))
        self.label_15 = QtGui.QLabel(Form)
        self.label_15.setEnabled(False)
        self.label_15.setGeometry(QtCore.QRect(870, 290, 181, 31))
        self.label_15.setAutoFillBackground(False)
        self.label_15.setWordWrap(True)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.lineEdit_16 = QtGui.QLineEdit(Form)
        self.lineEdit_16.setEnabled(False)
        self.lineEdit_16.setGeometry(QtCore.QRect(460, 480, 41, 31))
        self.lineEdit_16.setObjectName(_fromUtf8("lineEdit_16"))
        self.label_22 = QtGui.QLabel(Form)
        self.label_22.setEnabled(False)
        self.label_22.setGeometry(QtCore.QRect(30, 440, 121, 31))
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.lineEdit_19 = QtGui.QLineEdit(Form)
        self.lineEdit_19.setEnabled(False)
        self.lineEdit_19.setGeometry(QtCore.QRect(190, 480, 41, 31))
        self.lineEdit_19.setObjectName(_fromUtf8("lineEdit_19"))
        self.label_25 = QtGui.QLabel(Form)
        self.label_25.setEnabled(False)
        self.label_25.setGeometry(QtCore.QRect(30, 480, 181, 31))
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.label_26 = QtGui.QLabel(Form)
        self.label_26.setEnabled(False)
        self.label_26.setGeometry(QtCore.QRect(280, 480, 181, 31))
        self.label_26.setAutoFillBackground(False)
        self.label_26.setWordWrap(True)
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.lineEdit_11 = QtGui.QLineEdit(Form)
        self.lineEdit_11.setEnabled(False)
        self.lineEdit_11.setGeometry(QtCore.QRect(800, 380, 41, 31))
        self.lineEdit_11.setObjectName(_fromUtf8("lineEdit_11"))
        self.label_16 = QtGui.QLabel(Form)
        self.label_16.setEnabled(False)
        self.label_16.setGeometry(QtCore.QRect(30, 340, 111, 31))
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.label_17 = QtGui.QLabel(Form)
        self.label_17.setEnabled(False)
        self.label_17.setGeometry(QtCore.QRect(870, 380, 181, 31))
        self.label_17.setAutoFillBackground(False)
        self.label_17.setWordWrap(True)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.lineEdit_12 = QtGui.QLineEdit(Form)
        self.lineEdit_12.setEnabled(False)
        self.lineEdit_12.setGeometry(QtCore.QRect(1030, 380, 41, 31))
        self.lineEdit_12.setObjectName(_fromUtf8("lineEdit_12"))
        self.lineEdit_13 = QtGui.QLineEdit(Form)
        self.lineEdit_13.setEnabled(False)
        self.lineEdit_13.setGeometry(QtCore.QRect(150, 380, 41, 31))
        self.lineEdit_13.setObjectName(_fromUtf8("lineEdit_13"))
        self.lineEdit_14 = QtGui.QLineEdit(Form)
        self.lineEdit_14.setEnabled(False)
        self.lineEdit_14.setGeometry(QtCore.QRect(540, 380, 41, 31))
        self.lineEdit_14.setObjectName(_fromUtf8("lineEdit_14"))
        self.label_18 = QtGui.QLabel(Form)
        self.label_18.setEnabled(False)
        self.label_18.setGeometry(QtCore.QRect(220, 380, 71, 31))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.label_19 = QtGui.QLabel(Form)
        self.label_19.setEnabled(False)
        self.label_19.setGeometry(QtCore.QRect(370, 380, 181, 31))
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.lineEdit_15 = QtGui.QLineEdit(Form)
        self.lineEdit_15.setEnabled(False)
        self.lineEdit_15.setGeometry(QtCore.QRect(310, 380, 41, 31))
        self.lineEdit_15.setObjectName(_fromUtf8("lineEdit_15"))
        self.label_20 = QtGui.QLabel(Form)
        self.label_20.setEnabled(False)
        self.label_20.setGeometry(QtCore.QRect(610, 380, 181, 31))
        self.label_20.setAutoFillBackground(False)
        self.label_20.setWordWrap(True)
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.label_21 = QtGui.QLabel(Form)
        self.label_21.setEnabled(False)
        self.label_21.setGeometry(QtCore.QRect(30, 380, 111, 31))
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.label_6 = QtGui.QLabel(Form)
        self.label_6.setGeometry(QtCore.QRect(160, 90, 31, 31))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label_7 = QtGui.QLabel(Form)
        self.label_7.setGeometry(QtCore.QRect(410, 90, 31, 31))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_8 = QtGui.QLabel(Form)
        self.label_8.setGeometry(QtCore.QRect(610, 90, 31, 31))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_23 = QtGui.QLabel(Form)
        self.label_23.setEnabled(True)
        self.label_23.setGeometry(QtCore.QRect(20, 200, 1071, 31))
        self.label_23.setAutoFillBackground(False)
        self.label_23.setWordWrap(True)
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.label_24 = QtGui.QLabel(Form)
        self.label_24.setGeometry(QtCore.QRect(810, 550, 111, 31))
        self.label_24.setObjectName(_fromUtf8("label_24"))

        self.retranslateUi(Form)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_9.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_4.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_10.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_10.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_5.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_6.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_12.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_7.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_11.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.checkBox_4.setDisabled)
        QtCore.QObject.connect(self.checkBox_4, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.checkBox_2.setDisabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_16.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_13.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_21.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_13.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_18.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_15.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_19.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_14.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_20.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_11.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_17.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_12.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_15.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_10.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_22.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_25.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_19.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_26.setEnabled)
        QtCore.QObject.connect(self.checkBox_2, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_16.setEnabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_2.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_6.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_3.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_2.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_7.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_4.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.lineEdit_3.setDisabled)
        QtCore.QObject.connect(self.checkBox, QtCore.SIGNAL(_fromUtf8("toggled(bool)")), self.label_8.setDisabled)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        """
    This function is used to give names to the various components of the interface
        """
        Form.setWindowTitle(_translate("Form", "Lithium Ion Simulator: Input", None))
        self.label.setText(_translate("Form", "Select the width for Anode, Cathode and Separator ", None))
        self.label_2.setText(_translate("Form", "Anode", None))
        self.label_3.setText(_translate("Form", "Cathode", None))
        self.label_4.setText(_translate("Form", "Separator", None))
        self.checkBox.setText(_translate("Form", "Use Default Values", None))
        self.checkBox_2.setText(_translate("Form", "Input Your Own Parameters", None))
        self.label_9.setText(_translate("Form", "Conductivity", None))
        self.label_10.setText(_translate("Form", "Porosity", None))
        self.label_11.setText(_translate("Form", "Diffusion_Coeff", None))
        self.label_12.setText(_translate("Form", "Reaction Constant", None))
        self.checkBox_4.setText(_translate("Form", "Use Default Materials", None))
        self.label_5.setText(_translate("Form", "Duration of Simulation", None))
        self.checkBox_5.setText(_translate("Form", "Value of Current Density ", None))
        self.label_14.setText(_translate("Form", "sec", None))
        self.pushButton.setText(_translate("Form", "Solve", None))
        self.label_13.setText(_translate("Form", "For Anode", None))
        self.label_15.setText(_translate("Form", "Diff_Coeff_Solid", None))
        self.label_22.setText(_translate("Form", "For Separator", None))
        self.label_25.setText(_translate("Form", "Diffusion_Coeff", None))
        self.label_26.setText(_translate("Form", "Reaction Constant", None))
        self.label_16.setText(_translate("Form", "For Cathode", None))
        self.label_17.setText(_translate("Form", "Diff_Coeff_Solid", None))
        self.label_18.setText(_translate("Form", "Porosity", None))
        self.label_19.setText(_translate("Form", "Diffusion_Coeff", None))
        self.label_20.setText(_translate("Form", "Reaction Constant", None))
        self.label_21.setText(_translate("Form", "Conductivity", None))
        self.label_6.setText(_translate("Form", "m", None))
        self.label_7.setText(_translate("Form", "m", None))
        self.label_8.setText(_translate("Form", "m", None))
        self.label_23.setText(_translate("Form", "Standard Units: Conductivity = S/m, Porosity = None, Diffusion Coefficient = m^2/sec, Reaction Constant = (mol^5/2)/(s-(m^13/2))", None))
        self.label_24.setText(_translate("Form", "A/m^2", None))


    @QtCore.pyqtSignature("on_pushButton_clicked()")
    def Process(self):
	"""
	    	This function processes the input given by the user and assigns the input to the class attributes
	"""
        self.time = float(self.lineEdit_8.text())
        self.current_density = float(self.lineEdit_9.text())
        if not self.checkBox.isChecked():
            self.anode_width = float(self.lineEdit.text())
            self.cathode_width = float(self.lineEdit_2.text())
            self.separator_width = float(self.lineEdit_3.text())

        if self.checkBox_2.isChecked():

            #Taking input for anode material
            self.conductivity_anode = float(self.lineEdit_4.text())
            self.porosity_anode = float(self.lineEdit_5.text())
            self.duffusion_coeff_anode = float(self.lineEdit_6.text())
            self.duffusion_coeff__solid_anode = float(self.lineEdit_10.text())
            self.reaction_constant_anode = float(self.lineEdit_7.text())

            #Taking input for cathode material
            self.conductivity_cathode = float(self.lineEdit_13.text())
            self.porosity_cathode = float(self.lineEdit_15.text())
            self.duffusion_coeff_cathode = float(self.lineEdit_14.text())
            self.duffusion_coeff__solid_cathode = float(self.lineEdit_12.text())
            self.reaction_constant_cathode = float(self.lineEdit_11.text())

            #Taking input for separator
            self.reaction_constant_separator = float(self.lineEdit_16.text())
            self.duffusion_coeff_separator = float(self.lineEdit_25.text())




if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    ui = Ui_Form()
    ui.show()
    sys.exit(app.exec_())
