import json
import csv
from dataclasses import dataclass
from parse import *

@dataclass(frozen=True)
class EMRunTag:
    L: int
    N: int
    S: int
    E: float
    D: int
    method: str
    K: int

def parse_run_tag(s):
    parsed = parse('./data/L{L}_N{N}_S{S}/E{E}_D{D}/{method}_K{K}.tsv', s)
    return EMRunTag(
        L=int(parsed['L']),
        N=int(parsed['N']),
        S=int(parsed['S']),
        E=float(parsed['E']),
        D=int(parsed['D']),
        method=parsed['method'],
        K=int(parsed['K']),
    )

@dataclass
class EMLog:
    id: str
    depth: float
    prob: float
    seq: str
    copy_nums: list
    freqs: list

def parse_em_logs(tsv_filename):
    with open(tsv_filename) as f:
        logs = []
        for line in f:
            row = line.split('\t')
            log = EMLog(
                id=row[0],
                depth=float(row[1]),
                prob=float(row[2]),
                seq=row[3],
                copy_nums=json.loads(row[4]),
                freqs=json.loads(row[5]),
            )
            logs.append(log)
    return logs

@dataclass
class OptimizeLog:
    id: int
    temp: float
    now_score: float
    next_score: float
    p_accept: float
    is_accepted: bool
    now_score_prior: float
    now_score_forward: float
    now_size: int
    now_seq: str
    now_state: list
    next_score_prior: float
    next_score_forward: float
    next_size: int
    next_seq: str
    next_state: list

@dataclass
class GradComment:
    type: str
    score: float
    score_prior: float
    score_forward: float
    size: int
    seq: str
    state: list

def parse_logs(tsv_filename, with_comments=False):
    with open(tsv_filename) as f:
        logs = []
        comments = []
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0][0] == '#':
                if with_comments:
                    comment = GradComment(
                        type=row[0][2:],
                        score=float(row[1]),
                        score_prior=float(row[2]),
                        score_forward=float(row[3]),
                        size=int(row[4]),
                        seq=row[5],
                        state=json.loads(row[6]),
                    )
                    comments.append(comment)
            else:
                log = OptimizeLog(
                    id=int(row[0]),
                    temp=float(row[1]),
                    now_score=float(row[2]),
                    next_score=float(row[3]),
                    p_accept=float(row[4]),
                    is_accepted=(row[5] == 'true'),
                    now_score_prior=float(row[6]),
                    now_score_forward=float(row[7]),
                    now_size=int(row[8]),
                    now_seq=row[9],
                    now_state=json.loads(row[10]),
                    next_score_prior=float(row[11]),
                    next_score_forward=float(row[12]),
                    next_size=int(row[13]),
                    next_seq=row[14],
                    next_state=json.loads(row[15]),
                )
                logs.append(log)
    return logs, comments

def parse_stats(json_filename):
    with open(json_filename) as f:
        stats = json.load(f)
    return stats
