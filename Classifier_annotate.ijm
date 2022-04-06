//Notes:

//Annotate trash

//Debug:
//stop script here
//setBatchMode("exit and display");
//stop

//Variables:
start = 1; //starting point in the list

//Prepare workspace:
run("Close All");
print("\\Clear");

setBatchMode(true);

// Get paths
{
    dataDir = getDirectory("Process next files:");
    //print(dataDir)
    //dataDir = 'D:/Temp/APID/APID_5_Screen2/Morpho/8_SignatureAnalysis/CellImages/Output/Raw/';
    resDir = getDirectory("Save data to:");
    //resDir = 'D:/Temp/APID/APID_5_Screen2/Morpho/8_SignatureAnalysis/CellImages/Output/Raw/';
}
// Get file list
list = getFileList(dataDir);
//get cellN
cellN = lengthOf(list)
//array for decisions
annotation = newArray(cellN);
annotation_sure = newArray(cellN);

//Work through cells & annotate
for (cell = start-1; cell < cellN; cell+=1)
{
  //Path to img i
  img_path = dataDir+list[cell];

  //Get image Name
  FileNameFull = list[cell];
  dotIndex = lastIndexOf(FileNameFull, ".");
  if (dotIndex!=-1)
  FileName = substring(FileNameFull, 0, (dotIndex)); // name of file w/o extension


  //open image & move it to the stack
  open(img_path);
  //show to user
  setBatchMode("exit and display");
  run("In [+]");
  run("In [+]");
  run("In [+]");

  //Let user change how image looks
  waitForUser('Does cell look OK? If not -> adjust contrast etc.');

  sure=0;

  while (sure==0)
  {
    //Dialog with user
    Dialog.create("Annotator");
    Dialog.addMessage("Is this trash?");
    Dialog.addCheckbox("Yes", false);
    Dialog.addMessage("Are you sure?");
    Dialog.addCheckbox("Confirm choice", false);
    Dialog.show();

    //get answers
    answ1 = Dialog.getCheckbox();
    answ2 = Dialog.getCheckbox();
    //add answers to arrays
    annotation[cell] = answ1;
    annotation_sure[cell] = answ2;

    sure = answ2;
  }

  print(FileName+'/'+annotation[cell]);

  //Clean
  {
      //Closes all images:
      while (nImages>0) {
          selectImage(nImages);
          close();
      }
  }

}

//save result
selectWindow('Log');
save(resDir+'Result.txt');

//endpoint
//
//
//
//
//
//
//
//
