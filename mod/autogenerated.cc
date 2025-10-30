#include <cstdio>
#include <cstdlib>
#include <cstring>
void stntuple_print_header() {
  printf("stnmaker.exe built on XXDATE by $USER@hostname\n");
  printf("MU2E_BASE_RELEASE: release\n ");
  printf("OFFLINE  GIT tag : $offline_git_tag\n");
  printf("STNTUPLE GIT tag : $stntuple_git_tag\n");
}

void stntuple_get_version(char* Version, char* Text) {
  static char  txt[200];
  strcpy(Version,"$offline_version");
  strcpy(txt,"stnmaker.exe ");
  strcat(txt,"built on date by $USER@hostname");
  strcat(txt," MU2E_BASE_RELEASE $release, ");
  strcat(txt," OFFLINE  GIT tag $offline_git_tag  ");
  strcat(txt," STNTUPLE GIT tag $stntuple_git_tag ");
  strcpy(Text,txt);
}
