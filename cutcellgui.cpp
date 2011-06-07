/* Cartesian cut cell grid generator: GUI interface
 *
 * Given an OFF file and grid parameters, cut a grid and output the CGNS
 * representation of the resulting cartesian cut cell grid.
 *
 * Copyright 2010 Bruce Duncan, University of Edinburgh
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

    connect(cubicButton, SIGNAL(clicked()), this, SLOT(cubicGridSlot()));
    connect(offOpenButton, SIGNAL(clicked()), this, SLOT(getInFile()));
    connect(cgnsOpenButton, SIGNAL(clicked()), this, SLOT(getOutFile()));
    connect(generateButton, SIGNAL(clicked()), this, SLOT(generateSlot()));

    statusBarLabel = new QLabel();
    statusbar->addWidget(statusBarLabel);
}

void cutcellgui::cubicGridSlot() {
    double Xsize = gridSizeXSpinBox->value();
    gridSizeYSpinBox->setValue(Xsize);
    gridSizeZSpinBox->setValue(Xsize);
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
    thread = new GenerateThread(this);
    thread->start();
}

void GenerateThread::run() {
    cutcell::Polyhedron P;
    cutcellgui *p = reinterpret_cast<cutcellgui*>(parent());
    assert(p);
    // Store the values of the widgets' data now, before they can be changed.
    double transX = p->translationXSpinBox->value(),
           transY = p->translationYSpinBox->value(),
           transZ = p->translationZSpinBox->value(),
           scaling = p->scalingSpinBox->value();
    int gridX = p->gridSizeXSpinBox->value(),
        gridY = p->gridSizeYSpinBox->value(),
        gridZ = p->gridSizeZSpinBox->value();
    std::string outFile = qPrintable(p->cgnsFileLineEdit->text());
    // Read the OFF format from the input file or stdin.
    p->statusBarLabel->setText("Opening file...");
    std::ifstream in(qPrintable(p->offFileLineEdit->text()));
    if (!in) {
        // ifstream doesn't provide any way to retrieve the error message, though...
        p->statusBarLabel->setText("Open failed");
        return;
    }
    p->statusBarLabel->setText("Reading OFF...");
    // parse the file as a Nef_polyhedron object.
    CGAL::scan_OFF(in, P, true);
    if (in.bad()) {
        p->statusBarLabel->setText("Input file does not contain a permissible polyhedral surface.");
        return;
    }
    cutcell::Nef_polyhedron N1(P);
    // Transform the object accordingly.
    p->statusBarLabel->setText("Translating and scaling...");
    cutcell::Aff_transformation Aff1(cutcell::TRANSLATION,
                                     cutcell::Vector(transX, transY, transZ));
    cutcell::Aff_transformation Aff2(cutcell::SCALING, scaling);
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting grid.
    p->statusBarLabel->setText("Creating grid...");
    cutcell::Grid grid(gridX, gridY, gridZ);

    // Cut the grid to the object.
    p->statusBarLabel->setText("Cutting...");
    grid.addSolid(N1);
    grid.cut();

    // Output the grid in CGNS format.
    p->statusBarLabel->setText("Writing CGNS...");
    grid.output_cgns_file(outFile);
    p->statusBarLabel->setText("Done.");
}

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    cutcellgui *window = new cutcellgui;

    window->show();
    return app.exec();
}
