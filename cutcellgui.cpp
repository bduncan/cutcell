/* Cartesian cut cell grid generator: GUI interface
 *
 * Given an OFF file and grid parameters, cut a grid and output the CGNS
 * representation of the resulting cartesian cut cell grid.
 *
 * Copyright 2010,2011 Bruce Duncan, University of Edinburgh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cutcellgui.hpp>
#include <cutcell.hpp>
#include <QApplication>
#include <QtGui>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <string>

cutcellgui::cutcellgui(QWidget *parent) {
    setupUi(this);

    connect(cellSizeCubicButton, SIGNAL(clicked()), this, SLOT(cubicCellSlot()));
    connect(gridSizeCubicButton, SIGNAL(clicked()), this, SLOT(cubicGridSlot()));
    connect(originCubicButton, SIGNAL(clicked()), this, SLOT(cubicOriginSlot()));
    connect(offOpenButton, SIGNAL(clicked()), this, SLOT(getInFile()));
    connect(cgnsOpenButton, SIGNAL(clicked()), this, SLOT(getOutFile()));
    connect(generateButton, SIGNAL(clicked()), this, SLOT(generateSlot()));

    statusBarLabel = new QLabel();
    statusbar->addWidget(statusBarLabel);
}

void cutcellgui::cubicCellSlot() {
    double Xsize = cellSizeXDoubleSpinBox->value();
    cellSizeYDoubleSpinBox->setValue(Xsize);
    cellSizeZDoubleSpinBox->setValue(Xsize);
}

void cutcellgui::cubicGridSlot() {
    int Xsize = gridSizeXSpinBox->value();
    gridSizeYSpinBox->setValue(Xsize);
    gridSizeZSpinBox->setValue(Xsize);
}

void cutcellgui::cubicOriginSlot() {
    double Xsize = originXDoubleSpinBox->value();
    originYDoubleSpinBox->setValue(Xsize);
    originZDoubleSpinBox->setValue(Xsize);
}

void cutcellgui::getInFile() {
    QString path;

    path = QFileDialog::getOpenFileName(
        this,
        "Choose a file to open",
        QString::null,
        "OFF files (*.off)");

    offFileLineEdit->setText(path);
}

void cutcellgui::getOutFile() {
    QString path;

    path = QFileDialog::getSaveFileName(
        this,
        "Choose a file to write",
        QString::null,
        "CGNS files (*.cgns)");

    cgnsFileLineEdit->setText(path);
}

void cutcellgui::generateSlot() {
    // Disable the GUI so that the values cannot be changed as the thread reads them.
    setEnabled(false);
    thread = new GenerateThread(this);
    thread->start();
}

void GenerateThread::run() {
    cutcell::Polyhedron P;
    cutcellgui *p = reinterpret_cast<cutcellgui*>(parent());
    assert(p);
    // Read the OFF format from the input file.
    p->statusBarLabel->setText("Opening file...");
    std::ifstream in(qPrintable(p->offFileLineEdit->text()));
    if (!in) {
        // ifstream doesn't provide any way to retrieve the error message, though...
        p->statusBarLabel->setText("Open failed");
        p->setEnabled(true);
        return;
    }
    p->statusBarLabel->setText("Reading OFF...");
    // parse the file as a Nef_polyhedron object.
    CGAL::scan_OFF(in, P, true);
    if (in.bad()) {
        p->statusBarLabel->setText("Input file does not contain a permissible polyhedral surface.");
        p->setEnabled(true);
        return;
    }
    cutcell::Nef_polyhedron N1(P);
    // Transform the object accordingly.
    p->statusBarLabel->setText("Translating and scaling...");
    cutcell::Aff_transformation Aff1(cutcell::SCALING, p->scalingSpinBox->value());
    cutcell::Aff_transformation Aff2(cutcell::TRANSLATION,
                                     cutcell::Vector(p->translationXSpinBox->value(),
                                                     p->translationYSpinBox->value(),
                                                     p->translationZSpinBox->value()));
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting grid.
    p->statusBarLabel->setText("Creating grid...");
    cutcell::Grid grid(p->gridSizeXSpinBox->value(), p->gridSizeYSpinBox->value(), p->gridSizeZSpinBox->value());

    // Register the cell size and origin transform
    grid.addTransform(cutcell::Aff_transformation(p->cellSizeXDoubleSpinBox->value(), 0,  0,  p->originXDoubleSpinBox->value(),
                                                  0,  p->cellSizeYDoubleSpinBox->value(), 0,  p->originYDoubleSpinBox->value(),
                                                  0,  0,  p->cellSizeZDoubleSpinBox->value(), p->originZDoubleSpinBox->value()));
    // Cut the grid to the object.
    p->statusBarLabel->setText("Cutting...");
    grid.addSolid(N1);
    grid.cut();

    // Output the grid in CGNS format.
    p->statusBarLabel->setText("Writing CGNS...");
    grid.output_cgns_file(qPrintable(p->cgnsFileLineEdit->text()));
    p->statusBarLabel->setText("Done.");
    p->setEnabled(true);
}

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    cutcellgui *window = new cutcellgui;

    window->show();
    return app.exec();
}
