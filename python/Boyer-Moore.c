#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define BM_TABLE_SIZE 256

/*!
* @brief         ずらし表（照合再開位置テーブル）を作成する
* @param[in/out] table    テーブルへのアドレス
* @param[in]     pattern  探索文字列
* @param[in]     pat_len  パターン文字列長
*/
static void
bm_table_init(int *table, const char *pattern, int ptn_len)
{
	int cnt = 0;

	/* パターンに無い文字はパターン文字列長をずらし幅にする */
	for(cnt = 0; cnt < BM_TABLE_SIZE; cnt++){
		table[cnt] = ptn_len;
	}

	/* パターンに含まれる文字のずらし幅を決定する */
	for(cnt = 0; cnt < ptn_len; cnt++){
		table[(int)pattern[cnt]] = ptn_len - cnt - 1;
	}

	/* デバッグ出力 */
	printf("[table]  : default: step=%d\n", ptn_len);
	for(cnt = 0; cnt < BM_TABLE_SIZE; cnt++){
		if(table[cnt] != ptn_len)
			printf("         : char=%c: table[%03d]: step=%d\n",
					(char)cnt,cnt,(int)table[cnt]);
	}
}

/*!
* @brief    デバッグ出力
*/
static void
print_compare_process(const char *text, const char *pattern, int i, int j)
{
	int cnt = 0;

	printf("-----------------------------------\n");
	printf("[compare]:(text i=%d)(pattern j=%d)\n", i, j);
	printf(" text    :%s\n", text);

	/* パターンを比較位置で重ねる */
	printf(" pattern :");
	for(cnt = 0;cnt < (i - j); cnt++) printf(" ");
	printf("%s\n", pattern);

	/* 比較点にマークする */
	printf("         :");
	for(cnt = 0;cnt < i; cnt++) printf(" ");
	printf("^\n");
}

/*!
* @brief      パターンを進める加算値を求める
* @param[in]  table  ずらし表
* @param[in]  target 比較点の文字
* @param[in]  remain パターンの未比較部分の文字列長
* @return     パターンを進める加算値。
*/
static int
next_step(int *table, char target, int remain)
{
	/* ループ防止のために確認する */
	if(table[(int)target] > remain){
		/* ずらし表から値を取得する */
		return(table[(int)target]);
	}else{
		/* 照合を開始した地点の次の文字に進める */
		return(remain);
	}
}

/*!
* @brief      文字列を探索する
* @param[in]  text     検索対象文字列
* @param[in]  pattern  検索文字列
* @return     発見位置のポインタを返す。失敗したらNULLを返す。
*/
char *
bm_search(const char *text, const char *pattern)
{
	int table[BM_TABLE_SIZE];
	int txt_len = 0;
	int ptn_len = 0;
	int i = 0; /* テキストの比較位置 */
	int j = 0; /* パターンの比較位置 */

	ptn_len = strlen(pattern);
	txt_len = strlen(text);

	/* ずらし表を作成する */
	bm_table_init(table, pattern, ptn_len);

	/* 比較処理 */
	i = j = ptn_len - 1; /* 比較位置をパターン末尾にする */
	while((i < txt_len) && (j >= 0)){
		print_compare_process(text, pattern, i, j);

		if(text[i] != pattern[j]){
			/* ずらし表を参照して、次の比較位置を設定する */
			i += next_step(table, text[i], (ptn_len - j));
			j = ptn_len - 1;   /* 比較位置をパターン末尾にする */
		}else{
			/* 文字が一致したので、前の文字を照合していく */
			j--;
			i--;
		}
	}

	if(j < 0) return((char *)text + (i + 1));
	return(NULL);
}

int main(void) {
	char *text    = "GCTCACTGAGCGCTCGT";
	char *pattern = "GCTCG";
	char *cp = NULL;

	printf("[text]   :%s\n", text);
	printf("[pattern]:%s\n", pattern);

	cp = bm_search(text, pattern);
	if(cp == NULL){
		printf("[result] : not found\n");
	}else{
		printf("[result] : found\n");
		printf("         : text=%s\n", cp);
	}

	return(EXIT_SUCCESS);
}
