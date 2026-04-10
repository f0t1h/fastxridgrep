#ifndef FASTX_KSEQ_H
#define FASTX_KSEQ_H

#include <stddef.h>

typedef struct {
    void *gz;
    void *seq;
} FastxReader;

int fastx_open(FastxReader *reader, const char *path);
void fastx_close(FastxReader *reader);

int fastx_read(FastxReader *reader);

const char *fastx_name(const FastxReader *reader);
size_t fastx_name_len(const FastxReader *reader);

const char *fastx_comment(const FastxReader *reader);
size_t fastx_comment_len(const FastxReader *reader);

const char *fastx_seq(const FastxReader *reader);

const char *fastx_qual(const FastxReader *reader);
size_t fastx_qual_len(const FastxReader *reader);

#endif
