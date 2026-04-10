#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>

#include <zlib.h>

#include "../ext/aho-corasick/acism.h"
#include "../ext/kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    char **items;
    size_t *lens;
    MEMREF *refs;
    size_t count;
    size_t capacity;
} PatternStore;

typedef struct {
    const size_t *lens;
    size_t count;
    size_t id_len;
    bool matched;
} ScanContext;

static void print_usage(const char *prog)
{
    (void)fprintf(stderr,
                  "Usage: %s -p <id_file> [-i <fastx_file>]\n"
                  "\n"
                  "Options:\n"
                  "  -p, --patterns  Read IDs file (required)\n"
                  "  -i, --input     FASTA/FASTQ input path, '-' for stdin (default: -)\n"
                  "  -h, --help      Show this help message\n",
                  prog);
}

static void pattern_store_free(PatternStore *store)
{
    size_t i;

    if (store == NULL) {
        return;
    }

    for (i = 0U; i < store->count; i++) {
        free(store->items[i]);
    }

    free(store->items);
    free(store->lens);
    free(store->refs);

    store->items = NULL;
    store->lens = NULL;
    store->refs = NULL;
    store->count = 0U;
    store->capacity = 0U;
}

static int pattern_store_reserve(PatternStore *store, size_t need)
{
    size_t new_capacity;
    char **new_items;
    size_t *new_lens;

    if (store->capacity >= need) {
        return 0;
    }

    new_capacity = (store->capacity == 0U) ? 128U : store->capacity;
    while (new_capacity < need) {
        if (new_capacity > (SIZE_MAX / 2U)) {
            return -1;
        }
        new_capacity *= 2U;
    }

    new_items = (char **)realloc(store->items, new_capacity * sizeof(char *));
    if (new_items == NULL) {
        return -1;
    }

    new_lens = (size_t *)realloc(store->lens, new_capacity * sizeof(size_t));
    if (new_lens == NULL) {
        return -1;
    }

    store->items = new_items;
    store->lens = new_lens;
    store->capacity = new_capacity;
    return 0;
}

static int pattern_store_push(PatternStore *store, const char *token, size_t token_len)
{
    char *copy;

    if (pattern_store_reserve(store, store->count + 1U) != 0) {
        return -1;
    }

    copy = (char *)malloc(token_len + 1U);
    if (copy == NULL) {
        return -1;
    }

    if (token_len > 0U) {
        (void)memcpy(copy, token, token_len);
    }
    copy[token_len] = '\0';

    store->items[store->count] = copy;
    store->lens[store->count] = token_len;
    store->count++;

    return 0;
}

static int read_line_dynamic(FILE *fp, char **line_buf, size_t *line_len)
{
    size_t cap;
    size_t len;
    int ch;
    char *buf;
    char *tmp;

    cap = 256U;
    len = 0U;
    buf = (char *)malloc(cap);
    if (buf == NULL) {
        return -1;
    }

    for (;;) {
        ch = fgetc(fp);
        if (ch == EOF) {
            if (ferror(fp) != 0) {
                free(buf);
                return -1;
            }
            if (len == 0U) {
                free(buf);
                return 0;
            }
            break;
        }

        if ((char)ch == '\n') {
            break;
        }

        if ((len + 1U) >= cap) {
            size_t next_cap = cap * 2U;
            if (next_cap <= cap) {
                free(buf);
                return -1;
            }
            tmp = (char *)realloc(buf, next_cap);
            if (tmp == NULL) {
                free(buf);
                return -1;
            }
            buf = tmp;
            cap = next_cap;
        }

        buf[len] = (char)ch;
        len++;
    }

    buf[len] = '\0';
    *line_buf = buf;
    *line_len = len;
    return 1;
}

static int extract_first_token(const char *line, size_t line_len, const char **token, size_t *token_len)
{
    size_t i;

    if (line_len == 0U) {
        *token = NULL;
        *token_len = 0U;
        return 0;
    }

    if (isspace((unsigned char)line[0]) != 0) {
        *token = NULL;
        *token_len = 0U;
        return 0;
    }

    i = 0U;
    while ((i < line_len) && (isspace((unsigned char)line[i]) == 0)) {
        i++;
    }

    *token = line;
    *token_len = i;
    return 1;
}

static int load_patterns(const char *pattern_path, PatternStore *store)
{
    FILE *fp;
    int line_status;

    fp = fopen(pattern_path, "rb");
    if (fp == NULL) {
        (void)fprintf(stderr, "error: cannot open patterns file '%s': %s\n", pattern_path, strerror(errno));
        return -1;
    }

    for (;;) {
        char *line = NULL;
        size_t line_len = 0U;
        const char *token = NULL;
        size_t token_len = 0U;

        line_status = read_line_dynamic(fp, &line, &line_len);
        if (line_status < 0) {
            (void)fprintf(stderr, "error: failed reading patterns file '%s'\n", pattern_path);
            (void)fclose(fp);
            return -1;
        }

        if (line_status == 0) {
            break;
        }

        if ((line_len > 0U) && (line[line_len - 1U] == '\r')) {
            line[line_len - 1U] = '\0';
            line_len--;
        }

        if (extract_first_token(line, line_len, &token, &token_len) == 1) {
            if (pattern_store_push(store, token, token_len) != 0) {
                free(line);
                (void)fclose(fp);
                return -1;
            }
        }

        free(line);
    }

    if (fclose(fp) != 0) {
        (void)fprintf(stderr, "error: failed closing patterns file '%s'\n", pattern_path);
        return -1;
    }

    if (store->count == 0U) {
        (void)fprintf(stderr, "error: no valid read IDs found in '%s'\n", pattern_path);
        return -1;
    }

    store->refs = (MEMREF *)calloc(store->count, sizeof(MEMREF));
    if (store->refs == NULL) {
        (void)fprintf(stderr, "error: out of memory allocating pattern refs\n");
        return -1;
    }

    for (size_t i = 0U; i < store->count; i++) {
        store->refs[i].ptr = store->items[i];
        store->refs[i].len = store->lens[i];
    }

    return 0;
}

static int on_match_exact_id(int strnum, int textpos, void *context)
{
    ScanContext *ctx = (ScanContext *)context;
    size_t match_end;
    size_t pat_len;

    if ((ctx == NULL) || (strnum < 0) || (textpos < 0)) {
        return 0;
    }

    if ((size_t)strnum >= ctx->count) {
        return 0;
    }

    match_end = (size_t)textpos;
    pat_len = ctx->lens[(size_t)strnum];

    if ((pat_len == ctx->id_len) && (match_end == ctx->id_len)) {
        ctx->matched = true;
        return 1;
    }

    return 0;
}

static bool match_read_id_exact(ACISM *ac, const PatternStore *patterns, const char *id, size_t id_len)
{
    MEMREF text;
    ScanContext ctx;

    text.ptr = id;
    text.len = id_len;

    ctx.lens = patterns->lens;
    ctx.count = patterns->count;
    ctx.id_len = id_len;
    ctx.matched = false;

    (void)acism_scan(ac, text, on_match_exact_id, &ctx);

    return ctx.matched;
}

static int print_record(const kseq_t *rec, FILE *out)
{
    if (rec->qual.l > 0U) {
        if (rec->comment.l > 0U) {
            if (fprintf(out, "@%s %s\n%s\n+\n%s\n", rec->name.s, rec->comment.s, rec->seq.s, rec->qual.s) < 0) {
                return -1;
            }
        } else {
            if (fprintf(out, "@%s\n%s\n+\n%s\n", rec->name.s, rec->seq.s, rec->qual.s) < 0) {
                return -1;
            }
        }
    } else {
        if (rec->comment.l > 0U) {
            if (fprintf(out, ">%s %s\n%s\n", rec->name.s, rec->comment.s, rec->seq.s) < 0) {
                return -1;
            }
        } else {
            if (fprintf(out, ">%s\n%s\n", rec->name.s, rec->seq.s) < 0) {
                return -1;
            }
        }
    }

    return 0;
}

int main(int argc, char **argv)
{
    const char *pattern_path = NULL;
    const char *input_path = "-";
    PatternStore patterns = {0};
    ACISM *ac = NULL;
    gzFile in = NULL;
    kseq_t *seq = NULL;
    int ret;
    int exit_code = EXIT_FAILURE;

    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        }

        if (((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--patterns") == 0)) && ((i + 1) < argc)) {
            pattern_path = argv[++i];
            continue;
        }

        if (((strcmp(argv[i], "-i") == 0) || (strcmp(argv[i], "--input") == 0)) && ((i + 1) < argc)) {
            input_path = argv[++i];
            continue;
        }

        (void)fprintf(stderr, "error: unknown or incomplete argument '%s'\n", argv[i]);
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (pattern_path == NULL) {
        (void)fprintf(stderr, "error: -p/--patterns is required\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (load_patterns(pattern_path, &patterns) != 0) {
        goto cleanup;
    }

    if (patterns.count > (size_t)INT_MAX) {
        (void)fprintf(stderr, "error: too many patterns; max supported is %d\n", INT_MAX);
        goto cleanup;
    }

    ac = acism_create(patterns.refs, (int)patterns.count);
    if (ac == NULL) {
        (void)fprintf(stderr, "error: failed to build aho-corasick automaton\n");
        goto cleanup;
    }

    if (strcmp(input_path, "-") == 0) {
        in = gzdopen(fileno(stdin), "rb");
    } else {
        in = gzopen(input_path, "rb");
    }

    if (in == NULL) {
        (void)fprintf(stderr, "error: cannot open input '%s'\n", input_path);
        goto cleanup;
    }

    seq = kseq_init(in);
    if (seq == NULL) {
        (void)fprintf(stderr, "error: failed to initialize FASTX parser\n");
        goto cleanup;
    }

    for (;;) {
        ret = kseq_read(seq);
        if (ret < 0) {
            break;
        }

        if (match_read_id_exact(ac, &patterns, seq->name.s, seq->name.l)) {
            if (print_record(seq, stdout) != 0) {
                (void)fprintf(stderr, "error: failed writing output\n");
                goto cleanup;
            }
        }
    }

    if (ret != -1) {
        (void)fprintf(stderr, "error: malformed FASTX input (kseq code %d)\n", ret);
        goto cleanup;
    }

    if (fflush(stdout) != 0) {
        (void)fprintf(stderr, "error: failed flushing output\n");
        goto cleanup;
    }

    exit_code = EXIT_SUCCESS;

cleanup:
    if (seq != NULL) {
        kseq_destroy(seq);
    }
    if (in != NULL) {
        (void)gzclose(in);
    }
    if (ac != NULL) {
        acism_destroy(ac);
    }
    pattern_store_free(&patterns);

    return exit_code;
}
