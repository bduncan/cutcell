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
    // It's not clear that this provides actual thread-safety, but it seems to work in testing.
    setEnabled(false);
    thread = new GenerateThread(this);
    connect(thread, SIGNAL(done(bool, QString)), this, SLOT(done(bool, QString)));
    thread->start();
}

void cutcellgui::done(bool enabled, QString message) {
    statusBarLabel->setText(message);
    setEnabled(enabled);
}

void GenerateThread::run() {
    cutcell::Polyhedron P;
    cutcellgui *p = reinterpret_cast<cutcellgui*>(parent());
    assert(p);
    // A quick check that the output file has at least been set...
    if (p->cgnsFileLineEdit->text().isEmpty()) {
        emit done(true, "Output file not specified.");
        return;
    }
    // Read the OFF format from the input file.
    emit done(false, "Opening file...");
    std::ifstream in(qPrintable(p->offFileLineEdit->text()));
    if (!in) {
        // ifstream doesn't provide any way to retrieve the error message, though...
        emit done(true, "Open failed.");
        return;
    }
    emit done(false, "Reading OFF...");
    // parse the file as a Nef_polyhedron object.
    CGAL::scan_OFF(in, P, true);
    if (in.bad()) {
        emit done(true, "Input file does not contain a permissible polyhedral surface.");
        return;
    }
    cutcell::Nef_polyhedron N1(P);
    // Transform the object accordingly.
    emit done(false, "Translating and scaling...");
    cutcell::Aff_transformation Aff1(cutcell::SCALING, p->scalingSpinBox->value());
    cutcell::Aff_transformation Aff2(cutcell::TRANSLATION,
                                     cutcell::Vector(p->translationXSpinBox->value(),
                                                     p->translationYSpinBox->value(),
                                                     p->translationZSpinBox->value()));
    N1.transform(Aff1);
    N1.transform(Aff2);

    // Create the cutting grid.
    emit done(false, "Creating grid...");
    cutcell::Grid grid(p->gridSizeXSpinBox->value(), p->gridSizeYSpinBox->value(), p->gridSizeZSpinBox->value());

    // Register the cell size and origin transform
    grid.addTransform(cutcell::Aff_transformation(p->cellSizeXDoubleSpinBox->value(), 0,  0,  p->originXDoubleSpinBox->value(),
                                                  0,  p->cellSizeYDoubleSpinBox->value(), 0,  p->originYDoubleSpinBox->value(),
                                                  0,  0,  p->cellSizeZDoubleSpinBox->value(), p->originZDoubleSpinBox->value()));
    // Cut the grid to the object.
    emit done(false, "Cutting...");
    grid.addSolid(N1);
    grid.cut();

    // Output the grid in CGNS format.
    emit done(false, "Writing CGNS...");
    grid.output_cgns_file(qPrintable(p->cgnsFileLineEdit->text()));
    emit done(true, "Done.");
}

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    cutcellgui *window = new cutcellgui;

    window->show();
    return app.exec();
}
