use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

use clap::ValueEnum;
use indexmap::IndexMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum ClusterMethod {
    Unique,
    Percentile,
    Cluster,
    Adjacency,
    Directional,
}

pub fn cluster_umis(
    counts: &IndexMap<Vec<u8>, usize>,
    method: ClusterMethod,
    threshold: usize,
) -> Vec<Vec<Vec<u8>>> {
    let umis: Vec<Vec<u8>> = counts.keys().cloned().collect();
    let count_vec: Vec<usize> = umis.iter().map(|umi| counts[umi]).collect();

    match method {
        ClusterMethod::Unique => group_unique(&umis),
        ClusterMethod::Percentile => group_percentile(&umis, &count_vec),
        ClusterMethod::Cluster => {
            let adj = adjacency_graph(&umis, &count_vec, threshold, false);
            let components = connected_components(&umis, &count_vec, &adj);
            components
                .into_iter()
                .map(|mut component| {
                    sort_by_count(&mut component, &umis, &count_vec);
                    component.into_iter().map(|idx| umis[idx].clone()).collect()
                })
                .collect()
        }
        ClusterMethod::Adjacency => {
            let adj = adjacency_graph(&umis, &count_vec, threshold, false);
            let components = connected_components(&umis, &count_vec, &adj);
            group_adjacency(&umis, &count_vec, &adj, components)
        }
        ClusterMethod::Directional => {
            let adj = adjacency_graph(&umis, &count_vec, threshold, true);
            let components = connected_components(&umis, &count_vec, &adj);
            group_directional(&umis, &count_vec, components)
        }
    }
}

fn group_unique(umis: &[Vec<u8>]) -> Vec<Vec<Vec<u8>>> {
    if umis.len() == 1 {
        vec![umis.to_vec()]
    } else {
        umis.iter().cloned().map(|umi| vec![umi]).collect()
    }
}

fn group_percentile(umis: &[Vec<u8>], counts: &[usize]) -> Vec<Vec<Vec<u8>>> {
    if umis.len() == 1 {
        return vec![umis.to_vec()];
    }

    let threshold = median(counts) / 100.0;
    umis.iter()
        .zip(counts.iter())
        .filter(|(_, count)| **count as f64 > threshold)
        .map(|(umi, _)| vec![umi.clone()])
        .collect()
}

fn group_directional(
    umis: &[Vec<u8>],
    counts: &[usize],
    components: Vec<Vec<usize>>,
) -> Vec<Vec<Vec<u8>>> {
    let mut observed = HashSet::new();
    let mut groups = Vec::new();

    for mut component in components {
        if component.len() == 1 {
            observed.insert(component[0]);
            groups.push(vec![umis[component[0]].clone()]);
            continue;
        }

        sort_by_count(&mut component, umis, counts);
        let mut group = Vec::new();
        for idx in component {
            if observed.insert(idx) {
                group.push(umis[idx].clone());
            }
        }
        if !group.is_empty() {
            groups.push(group);
        }
    }

    groups
}

fn group_adjacency(
    umis: &[Vec<u8>],
    counts: &[usize],
    adj: &[Vec<usize>],
    components: Vec<Vec<usize>>,
) -> Vec<Vec<Vec<u8>>> {
    let mut groups = Vec::new();

    for component in components {
        if component.len() == 1 {
            groups.push(vec![umis[component[0]].clone()]);
            continue;
        }

        let lead_umis = best_min_account(&component, adj, counts, umis);
        let mut observed: HashSet<usize> = lead_umis.iter().copied().collect();

        for lead in lead_umis {
            let mut group = vec![umis[lead].clone()];
            let mut connected: Vec<usize> = adj[lead]
                .iter()
                .copied()
                .filter(|idx| !observed.contains(idx))
                .collect();
            sort_by_count(&mut connected, umis, counts);
            for idx in connected {
                observed.insert(idx);
                group.push(umis[idx].clone());
            }
            groups.push(group);
        }
    }

    groups
}

fn best_min_account(
    component: &[usize],
    adj: &[Vec<usize>],
    counts: &[usize],
    umis: &[Vec<u8>],
) -> Vec<usize> {
    if component.len() == 1 {
        return component.to_vec();
    }

    let mut sorted = component.to_vec();
    sort_by_count(&mut sorted, umis, counts);

    for end in 1..=sorted.len() {
        if remove_umis(adj, component, &sorted[..end]).is_empty() {
            return sorted[..end].to_vec();
        }
    }

    sorted
}

fn remove_umis(adj: &[Vec<usize>], component: &[usize], nodes: &[usize]) -> HashSet<usize> {
    let mut remove = HashSet::new();
    for node in nodes {
        remove.insert(*node);
        for neighbour in &adj[*node] {
            remove.insert(*neighbour);
        }
    }

    component
        .iter()
        .copied()
        .filter(|idx| !remove.contains(idx))
        .collect()
}

fn adjacency_graph(
    umis: &[Vec<u8>],
    counts: &[usize],
    threshold: usize,
    directional: bool,
) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); umis.len()];

    for (i, j) in candidate_pairs(umis, threshold) {
        if hamming_distance(&umis[i], &umis[j], threshold).is_some() {
            if directional {
                if counts[i] >= counts[j].saturating_mul(2).saturating_sub(1) {
                    adj[i].push(j);
                }
                if counts[j] >= counts[i].saturating_mul(2).saturating_sub(1) {
                    adj[j].push(i);
                }
            } else {
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }

    adj
}

fn connected_components(umis: &[Vec<u8>], counts: &[usize], adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut order: Vec<usize> = (0..umis.len()).collect();
    sort_by_count(&mut order, umis, counts);

    let mut found = HashSet::new();
    let mut components = Vec::new();

    for idx in order {
        if found.contains(&idx) {
            continue;
        }

        let mut component = breadth_first_search(idx, adj);
        component.sort_by(|a, b| umis[*a].cmp(&umis[*b]));
        found.extend(component.iter().copied());
        components.push(component);
    }

    components
}

fn breadth_first_search(node: usize, adj: &[Vec<usize>]) -> Vec<usize> {
    let mut searched = HashSet::new();
    let mut queue = VecDeque::new();
    searched.insert(node);
    queue.push_back(node);

    while let Some(current) = queue.pop_front() {
        for next in &adj[current] {
            if searched.insert(*next) {
                queue.push_back(*next);
            }
        }
    }

    searched.into_iter().collect()
}

fn candidate_pairs(umis: &[Vec<u8>], threshold: usize) -> Vec<(usize, usize)> {
    if umis.len() <= 25 {
        let mut pairs = Vec::new();
        for i in 0..umis.len() {
            for j in (i + 1)..umis.len() {
                pairs.push((i, j));
            }
        }
        return pairs;
    }

    let umi_len = umis.first().map_or(0, Vec::len);
    let slices = substring_slices(umi_len, threshold + 1);
    let mut index: Vec<HashMap<Vec<u8>, Vec<usize>>> =
        (0..slices.len()).map(|_| HashMap::new()).collect();

    for (slice_idx, (start, end)) in slices.iter().copied().enumerate() {
        for (umi_idx, umi) in umis.iter().enumerate() {
            index[slice_idx]
                .entry(umi[start..end].to_vec())
                .or_default()
                .push(umi_idx);
        }
    }

    let mut pairs = Vec::new();
    for (i, umi) in umis.iter().enumerate() {
        let mut neighbours = HashSet::new();
        for (slice_idx, (start, end)) in slices.iter().copied().enumerate() {
            if let Some(bucket) = index[slice_idx].get(&umi[start..end]) {
                for candidate in bucket {
                    if *candidate > i {
                        neighbours.insert(*candidate);
                    }
                }
            }
        }
        for j in neighbours {
            pairs.push((i, j));
        }
    }

    pairs
}

fn substring_slices(umi_length: usize, idx_size: usize) -> Vec<(usize, usize)> {
    if idx_size == 0 {
        return Vec::new();
    }

    let base = umi_length / idx_size;
    let remainder = umi_length % idx_size;
    let mut offset = 0;
    let mut slices = Vec::with_capacity(idx_size);

    for idx in 0..idx_size {
        let size = base + usize::from(idx < remainder);
        slices.push((offset, offset + size));
        offset += size;
    }

    slices
}

fn hamming_distance(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let mut distance = 0;
    for (x, y) in a.iter().zip(b.iter()) {
        if x != y {
            distance += 1;
            if distance > threshold {
                return None;
            }
        }
    }
    Some(distance)
}

fn sort_by_count(nodes: &mut [usize], umis: &[Vec<u8>], counts: &[usize]) {
    nodes.sort_by(|a, b| match counts[*b].cmp(&counts[*a]) {
        Ordering::Equal => umis[*a].cmp(&umis[*b]),
        ordering => ordering,
    });
}

fn median(values: &[usize]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }

    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    } else {
        sorted[mid] as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(items: &[(&[u8], usize)]) -> IndexMap<Vec<u8>, usize> {
        items
            .iter()
            .map(|(umi, count)| (umi.to_vec(), *count))
            .collect()
    }

    #[test]
    fn directional_groups_low_count_neighbour_under_parent() {
        let groups = cluster_umis(
            &counts(&[(b"AAAA", 10), (b"AAAT", 2), (b"TTTT", 4)]),
            ClusterMethod::Directional,
            1,
        );

        assert_eq!(groups[0], vec![b"AAAA".to_vec(), b"AAAT".to_vec()]);
        assert_eq!(groups[1], vec![b"TTTT".to_vec()]);
    }

    #[test]
    fn unique_keeps_each_umi_separate() {
        let groups = cluster_umis(
            &counts(&[(b"AAAA", 10), (b"AAAT", 2)]),
            ClusterMethod::Unique,
            1,
        );

        assert_eq!(groups, vec![vec![b"AAAA".to_vec()], vec![b"AAAT".to_vec()]]);
    }
}
