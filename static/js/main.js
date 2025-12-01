// 切換側邊欄顯示/隱藏
function toggleSidebar() {
    const sidebar = document.getElementById('sidebar');
    sidebar.classList.toggle('hidden');
}

function loadPage(pageName) {
    // 隱藏所有頁面
    document.querySelectorAll('.page-content').forEach(page => {
        page.classList.add('hidden');
    });

    // 顯示選中的頁面
    document.getElementById(pageName + '-page').classList.remove('hidden');

    // 更新導航按鈕樣式
    document.querySelectorAll('.nav-btn').forEach(btn => {
        btn.classList.remove('bg-blue-600', 'text-white');
        btn.classList.add('text-gray-300');
    });
    event.target.closest('.nav-btn').classList.add('bg-blue-600', 'text-white');
    event.target.closest('.nav-btn').classList.remove('text-gray-300');

    // 更新頁面標題 
    const titles = {
        'home': '首頁總覽',
        'normalize': '數據標準化',
        'umap': 'UMAP 視覺化',
        'heatmap': '基因表達熱圖',
        'annotation': '細胞類型註釋',
        'qc': 'QC 質控報告'
    };
    document.getElementById('page-title').textContent = titles[pageName];
}
async function apiRequest(endpoint, method = 'GET', data = null) {
    const options = {
        method: method,
        headers: {
            'Content-Type': 'application/json'
        }
    };
    
    if (data) {
        options.body = JSON.stringify(data);
    }
    
    try {
        const response = await fetch(endpoint, options);
        const result = await response.json();
        return result;
    } catch (error) {
        console.error('API request failed:', error);
        throw error;
    }
}

// 載入資料集列表
async function loadDatasets() {
    try {
        const datasets = await apiRequest('/api/datasets');
        return datasets;
    } catch (error) {
        console.error('Failed to load datasets:', error);
        return null;
    }
}

// 執行分析
async function runAnalysis(analysisType, params) {
    try {
        const result = await apiRequest('/api/analyze', 'POST', {
            type: analysisType,
            params: params
        });
        return result;
    } catch (error) {
        console.error('Analysis failed:', error);
        throw error;
    }
}

// 顯示載入中狀態
function showLoading(elementId, message = '處理中...') {
    const element = document.getElementById(elementId);
    if (element) {
        element.innerHTML = `
            <div class="flex items-center justify-center p-8">
                <svg class="animate-spin h-8 w-8 text-blue-500 mr-3" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                    <circle class="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" stroke-width="4"></circle>
                    <path class="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                </svg>
                <span class="text-gray-600">${message}</span>
            </div>
        `;
    }
}

// 顯示錯誤訊息
function showError(elementId, message) {
    const element = document.getElementById(elementId);
    if (element) {
        element.innerHTML = `
            <div class="bg-red-50 border border-red-200 rounded-lg p-4">
                <div class="flex items-center">
                    <svg class="w-5 h-5 text-red-500 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 8v4m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path>
                    </svg>
                    <span class="text-red-800">${message}</span>
                </div>
            </div>
        `;
    }
}

// 初始化
window.addEventListener('DOMContentLoaded', () => {
    // 選中首頁按鈕
    document.querySelector('.nav-btn').classList.add('bg-blue-600', 'text-white');
    
    // 載入資料集列表
    loadDatasets().then(datasets => {
        console.log('Datasets loaded:', datasets);
    });
});
