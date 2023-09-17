#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QString>
#include <QFileDialog>
#include <QFileInfo>

#include "rid-solomon.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , _input_zip_file(QFileInfo())
    , _input_unzip_file(QFileInfo())
    , _output_zip_file(QFileInfo())
    , _output_unzip_file(QFileInfo())
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_zip_clicked()
{
    FILE *input_file=nullptr, *output_file=nullptr;

    input_file=fopen(_input_zip_file.absoluteFilePath().toUtf8().constData(),"rb");
    output_file=fopen(_output_zip_file.absoluteFilePath().toUtf8().constData(),"wb");

    if (input_file==nullptr || output_file==nullptr){
        ui->listWidget_system_messages->addItem(QString("fopen: ошибка %1").arg(errno));
        ui->listWidget_system_messages->addItem(QString("Ошибка при открытии файла"));
    }else{
        int flag=code_file(input_file,_input_zip_file.absoluteFilePath().toUtf8().constData(),output_file);

        switch (flag) {
            case 0:
                ui->listWidget_system_messages->addItem(QString("Файл ")+_input_zip_file.absoluteFilePath()+ QString(" запакован в ")+_output_zip_file.absoluteFilePath());
            break;

            case 1:
                ui->listWidget_system_messages->addItem(QString("Ошибка при открытии файла"));
            break;

            case 2:
                ui->listWidget_system_messages->addItem(QString("Ошибка при выделении памяти"));
            break;

            case 3:
                ui->listWidget_system_messages->addItem(QString("Ошибка при перевыделении памяти"));
            break;

        }
    }
}


void MainWindow::on_pushButton_unzip_clicked()
{
    FILE *input_file=nullptr, *output_file=nullptr;

    input_file=fopen(_input_unzip_file.absoluteFilePath().toUtf8().constData(),"rb");
    output_file=fopen(_output_unzip_file.absoluteFilePath().toUtf8().constData(),"wb");

    if (input_file==nullptr || output_file==nullptr){
        ui->listWidget_system_messages->addItem(QString("fopen: ошибка %1").arg(errno));
        ui->listWidget_system_messages->addItem(QString("Ошибка при открытии файла"));
    }else{
        int flag=decode_file(input_file,output_file);

        switch (flag) {
            case 0:
                ui->listWidget_system_messages->addItem(QString("Файл ")+_input_unzip_file.absoluteFilePath()+ QString(" распакован в ")+_output_unzip_file.absoluteFilePath());
            break;

            case 1:
                ui->listWidget_system_messages->addItem(QString("Ошибка при открытии файла"));
            break;

            case 2:
                ui->listWidget_system_messages->addItem(QString("Ошибка при выделении памяти"));
            break;

            case 3:
                ui->listWidget_system_messages->addItem(QString("Ошибка при перевыделении памяти"));
            break;

        }
    }
}

void MainWindow::on_pushButton_input_zip_clicked()
{
    _input_zip_file.setFile(QFileDialog::getOpenFileName());
    if (_input_zip_file.absoluteFilePath()!=QString(""))
        ui->listWidget_system_messages->addItem(QString("Задан входной файл для запаковки ")+_input_zip_file.absoluteFilePath());
}

void MainWindow::on_pushButton_output_zip_clicked()
{
    _output_zip_file.setFile(QFileDialog::getSaveFileName());
    if (_output_zip_file.absoluteFilePath()!=QString(""))
        ui->listWidget_system_messages->addItem(QString("Задан выходной файл-архив для запаковки ")+_output_zip_file.absoluteFilePath());
}

void MainWindow::on_pushButton_input_unzip_clicked()
{
    _input_unzip_file.setFile(QFileDialog::getOpenFileName());
    if (_input_unzip_file.absoluteFilePath()!=QString(""))
        ui->listWidget_system_messages->addItem(QString("Задан входной файл-архив для распаковки ")+_input_unzip_file.absoluteFilePath());
}

void MainWindow::on_pushButton_output_unzip_clicked()
{
    _output_unzip_file.setFile(QFileDialog::getSaveFileName());
    if (_output_unzip_file.absoluteFilePath()!=QString(""))
        ui->listWidget_system_messages->addItem(QString("Задан выходной файл для распаковки ")+_output_unzip_file.absoluteFilePath());
}
