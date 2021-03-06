//===========================================================================
//===========================================================================
//===========================================================================
//==   Color2Gray_Sample.cpp  ==  Author: Diana-Maria Popa   ==
//===========================================================================
//===========================================================================
//===========================================================================

#include <sstream>
#include "stdafx.h"
#include "Color2Gray.h"
#include "ColorTransfer.h"
//===========================================================================
//===========================================================================
void runGOOCH();
void runPITIE();
void runFREEIMAGE();

inline const char * separator() {
#ifdef _WIN32
	return "\\";
#else
	return "/";
#endif
}

bool convertToChar(char *source, char **dest) {
	int size = 0;
	while ((char) source[size] != '\0') {
		size++;
	}
	size++;
	memcpy(*dest, source, size);
	return true;
}

enum Choice {
	GOOCH,
	PITIE,
	FREEIMAGE,
	NONE
};

bool removeExtension(char * str) {
	int len = strlen(str);

	len--;
	while (len >= 0 && str[len] != '.') {
		str[len] = 0;
		len--;
	}
	if (len >= 0) {
		str[len] = 0;
	}
	return true;
}

//basically you have to get the last part before the first slash

char * removeFolder(char * str) {
	char* token = strtok(str, separator());
	char* prevToken = token;
	while (token != NULL) {
		prevToken = token;
		token = strtok(NULL, separator());
	}
	return prevToken;
}

int main(int argc, char * argv[]) {
	if (argc < 2) {
		printf("You should provide at least the name of the input image!\n");
		return -1;
	}
	Choice processingType;
	switch (argc) {
		case 2:
			processingType = FREEIMAGE;
			printf("Applying FREEIMAGE...\n");
			break;
		case 3:
			processingType = PITIE;
			printf("Applying PITIE...\n");
			break;
		case 4:
			processingType = PITIE;
			printf("Applying PITIE...\n");
			break;
		case 5:
			processingType = GOOCH;
			printf("Applying GOOCH...\n");
			break;
		default:
			processingType = NONE;
			printf("Do not know what processing to apply !\n");
			return -1;
	}

	//Buffer for the new file names
	char strNewFileName[0x100];
	//Load and verify that input image is a True-Color one
	KImage *pImage = new KImage(argv[1]);
	if (pImage == NULL || !pImage->IsValid() || pImage->GetBPP() != 24) {
		printf("File %s does is not a valid True-Color image.\n", argv[1]);
		return -2;
	}
	if (processingType == GOOCH) {
		Color2Gray * clr2Gray = new Color2Gray(pImage);
		if (!clr2Gray) {
			printf("Could not allocate class Color2Gray. \n");
			return -1;
		}
		char buffer[1000];
		memset(buffer, 0, 1000);
		char * charpointer = new char[1000];

		convertToChar(argv[2], &charpointer);
		strcat(buffer, charpointer);
		strcat(buffer, " ");

		convertToChar(argv[3], &charpointer);
		strcat(buffer, charpointer);
		strcat(buffer, " ");

		convertToChar(argv[4], &charpointer);
		strcat(buffer, charpointer);

		//convert image with GOOCH algorithm

		KImage *pImageColor2Gray = clr2Gray->localC2G(buffer);
		sprintf(strNewFileName, "%s_GoochGray_%s_%s_%s.png", argv[1],
			argv[2], argv[3], argv[4]);
		pImageColor2Gray->SaveAs(strNewFileName, SAVE_PNG_DEFAULT);

		//convert image with BASIC algorithm
		KImage *pImageBasicGray = clr2Gray->localC2G(buffer,
			Color2Gray::BASIC_GREY);
		sprintf(strNewFileName, "%s_BasicGray.png", argv[1]);
		pImageBasicGray->SaveAs(strNewFileName, SAVE_PNG_DEFAULT);
	} else if (processingType == PITIE) {
		KImage *paletteImage = new KImage(argv[2]);
		if (paletteImage == NULL || !paletteImage->IsValid() ||
			paletteImage->GetBPP() != 24) {
			printf("Cannot perform Color Grading. \n");
			printf("File %s is not a valid True-Color image.",
				argv[2]);
			return -2;
		}
		ColorTransfer * Pitie = NULL;
		if (argc == 4) {
			char * pEnd = NULL;
			long int iterations = strtol(argv[3], &pEnd, 10);
			if (iterations == 0L) {
				printf("Invalid third argument for Pitie Color Grading.\n");
				printf("Provide a valid iterations number or accept the default value 10.\n");
				return -2;
			}
			Pitie = new ColorTransfer(pImage, paletteImage, iterations);
		} else {
			Pitie = new ColorTransfer(pImage, paletteImage);
		}
		assert(Pitie);
		KImage * recoloredImage = Pitie->applyRecoloring();
		assert(recoloredImage);
		removeExtension(argv[1]);
		argv[2] = removeFolder(argv[2]);
		removeExtension(argv[2]);
		sprintf(strNewFileName, "%s%s_Pitie.png", argv[1], argv[2]);
		printf("Result image saved to: %s\n", strNewFileName);
		recoloredImage->SaveAs(strNewFileName, SAVE_PNG_DEFAULT);
	} else if (processingType == FREEIMAGE) {
		//Convert to grayscale and then BlackAndWhite
		KImage *pImageGrayscale = pImage->ConvertToGrayscale();
		//Don't forget to delete the original, now useless image
		delete pImage;

		//Verify conversion success...
		if (pImageGrayscale == NULL || !pImageGrayscale->IsValid() || pImageGrayscale->GetBPP() != 8) {
			printf("Conversion to grayscale was not successful!");
			return -3;
		}
		//... and save grayscale image
		sprintf(strNewFileName, "%s_grayscale.TIF", argv[1]);
		pImageGrayscale->SaveAs(strNewFileName, SAVE_TIFF_LZW);
		//Request direct access to image pixels in raw format
		BYTE **pDataMatrixGrayscale = NULL;
		if (pImageGrayscale->BeginDirectAccess() &&
			(pDataMatrixGrayscale = pImageGrayscale->GetDataMatrix())
			!= NULL) {
			//If direct access is obtained get image attributes and
			//start processing pixels
			int intWidth = pImageGrayscale->GetWidth();
			int intHeight = pImageGrayscale->GetHeight();

			//Create binary image
			KImage *pImageBinary = new KImage(intWidth, intHeight, 1);
			if (pImageBinary->BeginDirectAccess()) {
				//Apply a threshold at half the grayscale range (0x00 is Full-Black, 0xFF is Full-White, 0x80 is the Middle-Gray)
				for (int y = intHeight - 1; y >= 0; y--)
					for (int x = intWidth - 1; x >= 0; x--) {
						//You may use this instead of the line below:
						//    BYTE PixelAtXY = pImageGrayscale->Get8BPPPixel(x, y)
						BYTE &PixelAtXY = pDataMatrixGrayscale[y][x];
						if (PixelAtXY < 0x80)
							//...if closer to black, set to black
							pImageBinary->Put1BPPPixel(x, y, false);
						else
							//...if closer to white, set to white
							pImageBinary->Put1BPPPixel(x, y, true);
					}

				//Close direct access
				pImageBinary->EndDirectAccess();

				//Save binarized image
				sprintf(strNewFileName, "%s_Black_and_White.TIF", argv[1]);
				pImageBinary->SaveAs(strNewFileName, SAVE_TIFF_CCITTFAX4);

				//Don't forget to delete the binary image
				delete pImageBinary;
			} else {
				printf("Unable to obtain direct access in binary image!");
				return -3;
			}

			//Close direct access
			pImageGrayscale->EndDirectAccess();
			delete pImageGrayscale;
		} else {
			printf("Unable to obtain direct access in grayscale image!");
			return -4;
		}
	}

	//Return with success
	return 0;
}
