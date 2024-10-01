#include <stdio.h>

void cti(float *x, float *y, int *n, FILE *fp);
float prumer(float *x, int n);

int main(){
    FILE *f;
    f = fopen("/home/reznicek/prednasky/comphep/cpp/xy.dat", "r");
    float x[20];
    float y[20];
    int n;

    cti(x, y, &n, f);

    float avarageX = prumer(x, n)/n;
    float avarageY = prumer(y, n)/n;

    printf("Avarage of x: %f\n", avarageX);
    printf("Avarage of y: %f\n", avarageY);

    fclose(f);
    return 0;
}

void cti(float *x, float *y, int *n, FILE *fp){
    *n = 0;

    while(fscanf(fp, "%f %f", &x[*n], &y[*n]) == 2){
        *n += 1;
    }
}

float prumer(float *x, int n){
    float sum = 0;
    for (int j = 0; j<n; j++){
        sum += x[j];
    }
    return sum;
}