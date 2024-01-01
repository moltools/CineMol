import React, { useState } from "react";

import Sidebar from "./Sidebar";

const WorkSpace = () => {
    const [sidebarOpen, setSideBarOpen] = useState(false);
    const handleViewSidebar = () => { setSideBarOpen(!sidebarOpen); };

    return (
        <div>
            <Sidebar isOpen={sidebarOpen} toggleSidebar={handleViewSidebar} />
        </div>
    );
};

export default WorkSpace;